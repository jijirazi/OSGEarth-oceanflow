/* -*-c++-*- */
/* osgEarth - Dynamic map generation toolkit for OpenSceneGraph
 * Copyright 2015 Pelican Mapping
 * http://osgearth.org
 *
 * osgEarth is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#include <osgEarthFeatures/BuildGeometryFilter>
#include <osgEarthFeatures/Session>
#include <osgEarthFeatures/FeatureSourceIndexNode>
#include <osgEarthFeatures/PolygonizeLines>
#include <osgEarthSymbology/TextSymbol>
#include <osgEarthSymbology/PointSymbol>
#include <osgEarthSymbology/LineSymbol>
#include <osgEarthSymbology/PolygonSymbol>
#include <osgEarthSymbology/MeshSubdivider>
#include <osgEarthSymbology/ResourceCache>
#include <osgEarth/Tessellator>
#include <osgEarth/Utils>
#include <osgEarth/Clamping>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/LineWidth>
#include <osg/LineStipple>
#include <osg/Point>
#include <osg/Depth>
#include <osg/PolygonOffset>
#include <osg/MatrixTransform>
#include <osgText/Text>
#include <osgUtil/Tessellator>
#include <osgUtil/Optimizer>
#include <osgUtil/Simplifier>
#include <osgUtil/SmoothingVisitor>
#include <osgDB/WriteFile>
#include <osg/Version>
#include <iterator>
#include <osgEarth/Registry>
#include <osgEarth/ShaderLoader>
#include <osgEarth/ObjectIndex>
#include <osgDB/WriteFile>
#include <osg/BlendFunc>
#include <osg/Texture1D>
#define LC "[BuildGeometryFilter] "

#define OE_TEST OE_NULL

using namespace osgEarth;
using namespace osgEarth::Features;
using namespace osgEarth::Symbology;

namespace
{
    bool isCCW(double x1, double y1, double x2, double y2, double x3, double y3)
    {
        return (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) > 0.0;
    }

    bool segmentsIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
    {
        return isCCW(x1, y1, x3, y3, x4, y4) != isCCW(x2, y2, x3, y3, x4, y4) && isCCW(x1, y1, x2, y2, x3, y3) != isCCW(x1, y1, x2, y2, x4, y4);
    }

    bool holeCompare(osgEarth::Symbology::Ring* i, osgEarth::Symbology::Ring* j)
    {
        return i->getBounds().xMax() > j->getBounds().xMax();
    }

    bool segmentsIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double &xi, double &yi)
    {
        double d = (y4-y3) * (x2-x1) - (x4-x3) * (y2-y1);

        if (d == 0) return false; // parallel

        double ua = ((x4-x3) * (y1-y3) - (y4-y3) * (x1-x3)) / d;
        double ub = ((x2-x1) * (y1-y3) - (y2-y1) * (x1-x3)) / d;

        if (ua >= 0.0 && ua <= 1.0 && ub >= 0.0 && ub <= 1.0)
        {
            xi = x1 + ua * (x2 - x1);
            yi = y1 + ua * (y2 - y1);

            return true;
        }

        return false;
    }
	const char* pickVertexEncode =
		"#version 130\n"

		"#pragma vp_entryPoint \"oe_pick_encodeObjectID\" \n"
		"#pragma vp_location   \"vertex_clip\" \n"

		"uint oe_index_objectid; \n"                        // Vertex stage global containing the Object ID; set in ObjectIndex shader.

		"flat out vec4 oe_pick_encoded_objectid; \n"        // output encoded oid to fragment shader
		"flat out int  oe_pick_color_contains_objectid; \n" // whether color already contains oid (written by another RTT camera)

		"void oe_pick_encodeObjectID(inout vec4 vertex) \n"
		"{ \n"
		"    oe_pick_color_contains_objectid = (oe_index_objectid == 1u) ? 1 : 0; \n"
		"    if ( oe_pick_color_contains_objectid == 0 ) \n"
		"    { \n"
		"        float b0 = float((oe_index_objectid & 0xff000000u) >> 24u); \n"
		"        float b1 = float((oe_index_objectid & 0x00ff0000u) >> 16u); \n"
		"        float b2 = float((oe_index_objectid & 0x0000ff00u) >> 8u ); \n"
		"        float b3 = float((oe_index_objectid & 0x000000ffu)       ); \n"
		"        oe_pick_encoded_objectid = vec4(b0, b1, b2, b3) * 0.00392156862; \n" // i.e. 1/2558
		"    } \n"
		"} \n";

	const char* pickFragment =
		"#version 130\n"

		"#pragma vp_entryPoint \"oe_pick_renderEncodedObjectID\" \n"
		"#pragma vp_location   \"fragment_output\" \n"
		"#pragma vp_order      \"last\" \n"

		"flat in vec4 oe_pick_encoded_objectid; \n"
		"flat in int  oe_pick_color_contains_objectid; \n"

		"out vec4 fragColor; \n"

		"void oe_pick_renderEncodedObjectID(inout vec4 color) \n"
		"{ \n"
		"    if ( oe_pick_color_contains_objectid == 1 ) \n"
		"        fragColor = color; \n"
		"    else \n"
		"        fragColor = oe_pick_encoded_objectid; \n"
		"} \n";
}
//
//float ZoomUpdateCallback::_scale = 1;

BuildGeometryFilter::BuildGeometryFilter( const Style& style ) :
_style        ( style ),
_maxAngle_deg ( 180.0 ),
_geoInterp    ( GEOINTERP_RHUMB_LINE )
{
    //nop
}

VirtualProgram* 
	BuildGeometryFilter::createRTTProgram()
{    
	VirtualProgram* vp = new VirtualProgram();
	vp->setName( "osgEarth::RTTPicker" );

	// Install RTT picker shaders:
	ShaderPackage pickShaders;
	pickShaders.add( "RTTPicker.vert.glsl", pickVertexEncode );
	pickShaders.add( "RTTPicker.frag.glsl", pickFragment );
	pickShaders.loadAll( vp );

	// Install shaders and bindings from the ObjectIndex:
	Registry::objectIndex()->loadShaders( vp );

	return vp;
}

osg::Geode*
BuildGeometryFilter::processPolygons(FeatureList& features, FilterContext& context)
{
    osg::Geode* geode = new osg::Geode();

    bool makeECEF = false;
    const SpatialReference* featureSRS = 0L;
    const SpatialReference* mapSRS = 0L;

    // set up the reference system info:
    if ( context.isGeoreferenced() )
    {
        makeECEF   = context.getSession()->getMapInfo().isGeocentric();
        featureSRS = context.extent()->getSRS();
        mapSRS     = context.getSession()->getMapInfo().getProfile()->getSRS();
    }

    for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
    {
        Feature* input = f->get();

        // access the polygon symbol, and bail out if there isn't one
        const PolygonSymbol* poly =
            input->style().isSet() && input->style()->has<PolygonSymbol>() ? input->style()->get<PolygonSymbol>() :
            _style.get<PolygonSymbol>();

        if ( !poly )
            continue;

        // run a symbol script if present.
        if ( poly->script().isSet() )
        {
            StringExpression temp( poly->script().get() );
            input->eval( temp, &context );
        }

        GeometryIterator parts( input->getGeometry(), false );
        while( parts.hasMore() )
        {
            Geometry* part = parts.next();

            part->removeDuplicates();

            // skip geometry that is invalid for a polygon
            if ( part->size() < 3 )
                continue;

            // resolve the color:
            osg::Vec4f primaryColor = poly->fill()->color();
            
            osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
            //osgGeom->setUseVertexBufferObjects( true );
            //osgGeom->setUseDisplayList( false );

            // are we embedding a feature name?
            if ( _featureNameExpr.isSet() )
            {
                const std::string& name = input->eval( _featureNameExpr.mutable_value(), &context );
                osgGeom->setName( name );
            }


            // compute localizing matrices or use globals
            osg::Matrixd w2l, l2w;
            if (makeECEF)
            {
                osgEarth::GeoExtent featureExtent(featureSRS);
                featureExtent.expandToInclude(part->getBounds());

                computeLocalizers(context, featureExtent, w2l, l2w);
            }
            else
            {
                w2l = _world2local;
                l2w = _local2world;
            }


            // build the geometry:
            tileAndBuildPolygon(part, featureSRS, mapSRS, makeECEF, true, osgGeom, w2l);
            //buildPolygon(part, featureSRS, mapSRS, makeECEF, true, osgGeom, w2l);

            osg::Vec3Array* allPoints = static_cast<osg::Vec3Array*>(osgGeom->getVertexArray());
            if (allPoints && allPoints->size() > 0)
            {
                // subdivide the mesh if necessary to conform to an ECEF globe:
                if ( makeECEF )
                {
                    //convert back to world coords
                    for( osg::Vec3Array::iterator i = allPoints->begin(); i != allPoints->end(); ++i )
                    {
                        osg::Vec3d v(*i);
                        v = v * l2w;
                        v = v * _world2local;

                        (*i)._v[0] = v[0];
                        (*i)._v[1] = v[1];
                        (*i)._v[2] = v[2];
                    }

                    double threshold = osg::DegreesToRadians( *_maxAngle_deg );
                    OE_DEBUG << "Running mesh subdivider with threshold " << *_maxAngle_deg << std::endl;

                    MeshSubdivider ms( _world2local, _local2world );
                    if ( input->geoInterp().isSet() )
                        ms.run( *osgGeom, threshold, *input->geoInterp() );
                    else
                        ms.run( *osgGeom, threshold, *_geoInterp );
                }

                // assign the primary color array. PER_VERTEX required in order to support
                // vertex optimization later
                unsigned count = osgGeom->getVertexArray()->getNumElements();
                osg::Vec4Array* colors = new osg::Vec4Array;
                colors->assign( count, primaryColor );
                osgGeom->setColorArray( colors );
                osgGeom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );

                geode->addDrawable( osgGeom );

                // record the geometry's primitive set(s) in the index:
                if ( context.featureIndex() )
                    context.featureIndex()->tagDrawable( osgGeom, input );
        
                // install clamping attributes if necessary
                if (_style.has<AltitudeSymbol>() &&
                    _style.get<AltitudeSymbol>()->technique() == AltitudeSymbol::TECHNIQUE_GPU)
                {            
                    Clamping::applyDefaultClampingAttrs( osgGeom, input->getDouble("__oe_verticalOffset", 0.0) );
                }
            }
        }
    }
    
    return geode;
}


osg::Geode*
BuildGeometryFilter::processPolygonizedLines(FeatureList&   features, 
                                             bool           twosided,
                                             FilterContext& context)
{
    osg::Geode* geode = new osg::Geode();

    // establish some referencing
    bool                    makeECEF   = false;
    const SpatialReference* featureSRS = 0L;
    const SpatialReference* mapSRS     = 0L;

    if ( context.isGeoreferenced() )
    {
        makeECEF   = context.getSession()->getMapInfo().isGeocentric();
        featureSRS = context.extent()->getSRS();
        mapSRS     = context.getSession()->getMapInfo().getProfile()->getSRS();
    }

    // iterate over all features.
    for( FeatureList::iterator i = features.begin(); i != features.end(); ++i )
    {
        Feature* input = i->get();
        // extract the required line symbol; bail out if not found.
        const LineSymbol* line =
            input->style().isSet() && input->style()->has<LineSymbol>() ? input->style()->get<LineSymbol>() :
            _style.get<LineSymbol>();

        if ( !line )
            continue;

        // run a symbol script if present.
        if ( line->script().isSet() )
        {
            StringExpression temp( line->script().get() );
            input->eval( temp, &context );
        }

        // The operator we'll use to make lines into polygons.
        PolygonizeLinesOperator polygonizer( *line->stroke() );

        // iterate over all the feature's geometry parts. We will treat
        // them as lines strings.
        GeometryIterator parts( input->getGeometry(), true );
        while( parts.hasMore() )
        {
            Geometry* part = parts.next();

            // if the underlying geometry is a ring (or a polygon), close it so the
            // polygonizer will generate a closed loop.
            Ring* ring = dynamic_cast<Ring*>(part);
            if ( ring )
                ring->close();

            // skip invalid geometry
            if ( part->size() < 2 )
                continue;

            // transform the geometry into the target SRS and localize it about 
            // a local reference point.
            osg::ref_ptr<osg::Vec3Array> verts   = new osg::Vec3Array();
            osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array();
            transformAndLocalize( part->asVector(), featureSRS, verts.get(), normals.get(), mapSRS, _world2local, makeECEF );

            // turn the lines into polygons.
            osg::Geometry* geom = polygonizer( verts.get(), normals.get(), twosided );
            if ( geom )
            {
                geode->addDrawable( geom );
            }

            // record the geometry's primitive set(s) in the index:
            if ( context.featureIndex() )
                context.featureIndex()->tagDrawable( geom, input );
        
            // install clamping attributes if necessary
            if (_style.has<AltitudeSymbol>() &&
                _style.get<AltitudeSymbol>()->technique() == AltitudeSymbol::TECHNIQUE_GPU)
            {
                Clamping::applyDefaultClampingAttrs( geom, input->getDouble("__oe_verticalOffset", 0.0) );
            }
        }

        polygonizer.installShaders( geode );
    }
    return geode;
}


osg::Geode*
BuildGeometryFilter::processLines(FeatureList& features, FilterContext& context)
{
    osg::Geode* geode = new osg::Geode();

    bool makeECEF = false;
    const SpatialReference* featureSRS = 0L;
    const SpatialReference* mapSRS = 0L;

    // set up referencing information:
    if ( context.isGeoreferenced() )
    {
        makeECEF   = context.getSession()->getMapInfo().isGeocentric();
        featureSRS = context.extent()->getSRS();
        mapSRS     = context.getSession()->getMapInfo().getProfile()->getSRS();
    }

    for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
    {
        Feature* input = f->get();

        // extract the required line symbol; bail out if not found.
        const LineSymbol* line = 
            input->style().isSet() && input->style()->has<LineSymbol>() ? input->style()->get<LineSymbol>() :
            _style.get<LineSymbol>();

        if ( !line )
            continue;

        // run a symbol script if present.
        if ( line->script().isSet() )
        {
            StringExpression temp( line->script().get() );
            input->eval( temp, &context );
        }

        GeometryIterator parts( input->getGeometry(), true );
        while( parts.hasMore() )
        {
            Geometry* part = parts.next();

            // skip invalid geometry for lines.
            if ( part->size() < 2 )
                continue;

            // if the underlying geometry is a ring (or a polygon), use a line loop; otherwise
            // use a line strip.
            GLenum primMode = dynamic_cast<Ring*>(part) ? GL_LINE_LOOP : GL_LINE_STRIP;

            // resolve the color:
            osg::Vec4f primaryColor = line->stroke()->color();
            
            osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
            //osgGeom->setUseVertexBufferObjects( true );
            //osgGeom->setUseDisplayList( false );

            // embed the feature name if requested. Warning: blocks geometry merge optimization!
            if ( _featureNameExpr.isSet() )
            {
                const std::string& name = input->eval( _featureNameExpr.mutable_value(), &context );
                osgGeom->setName( name );
            }

            // build the geometry:
            osg::Vec3Array* allPoints = new osg::Vec3Array();

            transformAndLocalize( part->asVector(), featureSRS, allPoints, mapSRS, _world2local, makeECEF );

            osgGeom->addPrimitiveSet( new osg::DrawArrays(primMode, 0, allPoints->getNumElements()) );
            osgGeom->setVertexArray( allPoints );

            if ( input->style().isSet() )
            {
                //TODO: re-evaluate this. does it hinder geometry merging?
                applyLineSymbology( osgGeom->getOrCreateStateSet(), line );
            }
            
            // subdivide the mesh if necessary to conform to an ECEF globe;
            // but if the tessellation is set to zero, or if the style specifies a
            // tessellation size, skip this step.
            if ( makeECEF && !line->tessellation().isSetTo(0) && !line->tessellationSize().isSet() )
            {
                double threshold = osg::DegreesToRadians( *_maxAngle_deg );
                OE_DEBUG << "Running mesh subdivider with threshold " << *_maxAngle_deg << std::endl;

                MeshSubdivider ms( _world2local, _local2world );
                //ms.setMaxElementsPerEBO( INT_MAX );
                if ( input->geoInterp().isSet() )
                    ms.run( *osgGeom, threshold, *input->geoInterp() );
                else
                    ms.run( *osgGeom, threshold, *_geoInterp );
            }

            // assign the primary color (PER_VERTEX required for later optimization)
            osg::Vec4Array* colors = new osg::Vec4Array;
            colors->assign( osgGeom->getVertexArray()->getNumElements(), primaryColor );
            osgGeom->setColorArray( colors );
            osgGeom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );

            geode->addDrawable( osgGeom );

            // record the geometry's primitive set(s) in the index:
            if ( context.featureIndex() )
                context.featureIndex()->tagDrawable( osgGeom, input );
        
            // install clamping attributes if necessary
            if (_style.has<AltitudeSymbol>() &&
                _style.get<AltitudeSymbol>()->technique() == AltitudeSymbol::TECHNIQUE_GPU)
            {
                Clamping::applyDefaultClampingAttrs( osgGeom, input->getDouble("__oe_verticalOffset", 0.0) );
            }
        }
    }
    
    return geode;
	//osg::Geode* geode = new osg::Geode();

 //   bool makeECEF = false;
 //   const SpatialReference* featureSRS = 0L;
 //   const SpatialReference* mapSRS = 0L;
	//
 //   // set up referencing information:
 //   if ( context.isGeoreferenced() )
 //   {
 //       makeECEF   = context.getSession()->getMapInfo().isGeocentric();
 //       featureSRS = context.extent()->getSRS();
 //       mapSRS     = context.getSession()->getMapInfo().getProfile()->getSRS();
 //   }
	////float nowState = 1.0;

	//int colorStep = 0;
	//
	//if(_style.get<LineSymbol>()->trackLineString().get().size()>0)
	//{
	//	osg::Vec3Array* allPoints = new osg::Vec3Array();	
	//	const LineSymbol* line = _style.get<LineSymbol>();
	//	osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
	//	osgGeom->setUseVertexBufferObjects( true );
	//	osgGeom->setUseDisplayList( false );
	//	osg::Vec4f primaryColor = line->stroke()->color();
	//	applyLineSymbology( osgGeom->getOrCreateStateSet(), line );
	//	int offset =0;
	//	//bool needDoSubdiv = false;
	//	


	//	int colorStep = 0;
	//	osg::Vec4Array* colorArray = new osg::Vec4Array;
	//	for(int k = 0;k<(*features.begin()).get()->getGeometry()->size();k++)
	//	{
	//		int colorStart = line->trackLineString().get().find('*',colorStep);
	//		int colorEnd = line->trackLineString().get().find('*',colorStart+1);
	//		colorStep = colorEnd;

	//		float trueColorFloat[3];
	//		std::string tempVecString = line->trackLineString().get().substr(colorStart+1,colorEnd-colorStart-1);
	//		{
	//			int tempStep = 0;

	//			for(int i=0;i<3;i++)
	//			{
	//				int tempStart = tempVecString.find(',',tempStep);
	//				int tempEnd = tempVecString.find(',',tempStart+1);
	//				tempStep = tempEnd;
	//				std::string tempColorCell = tempVecString.substr(tempStart+1,tempEnd-tempStart-1);
	//				std::istringstream iss(tempColorCell);
	//				iss>>trueColorFloat[i];
	//			}
	//		}

	//		osg::Vec4 trueColor = osg::Vec4(trueColorFloat[0],trueColorFloat[1],trueColorFloat[2],1.0);
	//		colorArray->push_back(trueColor);
	//	}


	//	std::vector<osg::Vec3d> vertexInput;
	//	for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
	//	{
	//		Feature* input = f->get();

	//		GeometryIterator parts( input->getGeometry(), true );
	//		while( parts.hasMore() )
	//		{
	//			Geometry* part = parts.next();
	//			if ( part->size() < 2 )
	//				continue;
	//			GLenum primMode = dynamic_cast<Ring*>(part) ? GL_LINE_LOOP : GL_LINE_STRIP;
	//			vertexInput.insert(vertexInput.end(),part->asVector().begin(), part->asVector().end());
	//			osgGeom->addPrimitiveSet( new osg::DrawArrays(GL_POINTS, offset, part->asVector().size()) );
	//			offset+=part->asVector().size();
	//		}
	//	}
	//	transformAndLocalize(vertexInput, featureSRS, allPoints, mapSRS, _world2local, makeECEF );

	//	osgGeom->setVertexArray( allPoints );
	//	//if ( makeECEF && !line->tessellation().isSetTo(0)/*&&needDoSubdiv*/)
	//	//{
	//	//	double threshold = osg::DegreesToRadians( *_maxAngle_deg );
	//	//	OE_DEBUG << "Running mesh subdivider with threshold " << *_maxAngle_deg << std::endl;

	//	//	MeshSubdivider ms( _world2local, _local2world );

	//	//	ms.run( *osgGeom, threshold, *_geoInterp );
	//	//}
	//	////osg::Vec4Array* colors = new osg::Vec4Array;
	//	////colors->assign( osgGeom->getVertexArray()->getNumElements(), primaryColor );
	//	//osgGeom->setColorArray( colorArray );
	//	//osgGeom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );

	//	//着色器 jxf
	//	osg::StateSet *ss = osgGeom->getOrCreateStateSet();
	//	osg::Program *program = new osg::Program;
	//	program->setParameter( GL_GEOMETRY_VERTICES_OUT_EXT, 4 );
	//	program->setParameter( GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS );
	//	program->setParameter( GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_LINE_STRIP );
	//	char* fragShader;
	//	char* vertShader;
	//	char* geoShader;

	//	osg::Uniform* u_anim1( new osg::Uniform( "u_anim1", 0.0f ) );
	//	u_anim1->setUpdateCallback( new SineAnimation( 4, 0.5, 0.5 ) );
	//	ss->addUniform( u_anim1 );

	//	vertShader = 
	//		"#version 330 compatibility\n"
	//		"uniform float u_anim1;\n"
	//		"varying vec4 v_color;\n"
	//		"void main(void)\n"
	//		"{\n" 
	//		"    v_color = vec4(1.0,0.0,0.0, 1.0);\n"
	//		"    gl_Position = gl_Vertex;\n"
	//		"}\n"
	//	;

	//	geoShader = 
	//		"#version 330 compatibility\n"
	//		"#extension GL_EXT_geometry_shader4 : enable\n"
	//		"uniform float u_anim1;\n"
	//		"layout(points) in;\n"
	//		"layout(line_strip,max_vertices=4) out;\n"
	//		"varying in vec4 v_color[];\n"
	//		"varying out vec4 v_color_out;\n"
	//		"void main(void)\n"
	//		"{\n"
	//		"    vec4 v = gl_PositionIn[0];\n"
	//		"    v_color_out = v_color[0];\n"
	//		"\n"
	//		/*"    gl_Position = gl_ModelViewProjectionMatrix * (v + vec4(u_anim1,0.,0.,0.));  EmitVertex();\n"
	//		"    gl_Position = gl_ModelViewProjectionMatrix * (v - vec4(u_anim1,0.,0.,0.));  EmitVertex();\n"*/
	//		"    gl_Position = gl_ModelViewProjectionMatrix * v;  EmitVertex();\n"
	//		"    gl_Position = gl_ModelViewProjectionMatrix * (v+vec4(10.0,0.,0.,0.));  EmitVertex();\n"
	//		/*"    EndPrimitive();\n"
	//		"\n"*/
	//		/*"    gl_Position = gl_ModelViewProjectionMatrix * (v + vec4(0.,1.0-u_anim1,0.,0.));  EmitVertex();\n"
	//		"    gl_Position = gl_ModelViewProjectionMatrix * (v - vec4(0.,1.0-u_anim1,0.,0.));  EmitVertex();\n"*/
	//		"    gl_Position = gl_ModelViewProjectionMatrix * (v+vec4(10.0,10.0,0.,0.));  EmitVertex();\n"
	//		"    gl_Position = gl_ModelViewProjectionMatrix * (v+vec4(0.,10.0,0.,0.));  EmitVertex();\n"
	//		"    EndPrimitive();\n"
	//		"}\n"
	//	;


	//	fragShader = 
	//		"#version 330 compatibility\n"
	//		"uniform float u_anim1;\n"
	//		"varying vec4 v_color_out;\n"
	//		"void main(void)\n"
	//		"{\n"
	//		"    gl_FragColor = v_color_out;\n"
	//		"}\n"                   
	//	;
	//	program->addShader(new osg::Shader(osg::Shader::FRAGMENT,fragShader));
	//	program->addShader(new osg::Shader(osg::Shader::VERTEX,vertShader));
	//	program->addShader(new osg::Shader(osg::Shader::GEOMETRY,geoShader));
	//	ss->setAttributeAndModes(program,osg::StateAttribute::ON);

	//	geode->addDrawable( osgGeom );
	//}
	//else
	//{
	//	for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
	//	{
	//	    Feature* input = f->get();

	//	    // extract the required line symbol; bail out if not found.
	//	    const LineSymbol* line = 
	//	        input->style().isSet() && input->style()->has<LineSymbol>() ? input->style()->get<LineSymbol>() :
	//	        _style.get<LineSymbol>();

	//	    if ( !line )
	//	        continue;

	//	    // run a symbol script if present.
	//	    if ( line->script().isSet() )
	//	    {
	//	        StringExpression temp( line->script().get() );
	//	        input->eval( temp, &context );
	//	    }

	//	    GeometryIterator parts( input->getGeometry(), true );
	//	    while( parts.hasMore() )
	//	    {
	//	        Geometry* part = parts.next();

	//	        // skip invalid geometry for lines.
	//	        if ( part->size() < 2 )
	//	            continue;

	//	        // if the underlying geometry is a ring (or a polygon), use a line loop; otherwise
	//	        // use a line strip.
	//	        GLenum primMode = dynamic_cast<Ring*>(part) ? GL_LINE_LOOP : GL_LINE_STRIP;

	//	        // resolve the color:
	//	        osg::Vec4f primaryColor = line->stroke()->color();
	//	        
	//	        osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
	//	        //osgGeom->setUseVertexBufferObjects( true );
	//	        //osgGeom->setUseDisplayList( false );

	//	        // embed the feature name if requested. Warning: blocks geometry merge optimization!
	//	        if ( _featureNameExpr.isSet() )
	//	        {
	//	            const std::string& name = input->eval( _featureNameExpr.mutable_value(), &context );
	//	            osgGeom->setName( name );
	//	        }

	//	        // build the geometry:
	//	        osg::Vec3Array* allPoints = new osg::Vec3Array();

	//	        transformAndLocalize( part->asVector(), featureSRS, allPoints, mapSRS, _world2local, makeECEF );

	//	        osgGeom->addPrimitiveSet( new osg::DrawArrays(primMode, 0, allPoints->getNumElements()) );
	//	        osgGeom->setVertexArray( allPoints );

	//	        if ( input->style().isSet() )
	//	        {
	//	            //TODO: re-evaluate this. does it hinder geometry merging?
	//	            applyLineSymbology( osgGeom->getOrCreateStateSet(), line );
	//	        }
	//	        
	//	        // subdivide the mesh if necessary to conform to an ECEF globe;
	//	        // but if the tessellation is set to zero, or if the style specifies a
	//	        // tessellation size, skip this step.
	//	        if ( makeECEF && !line->tessellation().isSetTo(0) && !line->tessellationSize().isSet() )
	//	        {
	//	            double threshold = osg::DegreesToRadians( *_maxAngle_deg );
	//	            OE_DEBUG << "Running mesh subdivider with threshold " << *_maxAngle_deg << std::endl;

	//	            MeshSubdivider ms( _world2local, _local2world );
	//	            //ms.setMaxElementsPerEBO( INT_MAX );
	//	            if ( input->geoInterp().isSet() )
	//	                ms.run( *osgGeom, threshold, *input->geoInterp() );
	//	            else
	//	                ms.run( *osgGeom, threshold, *_geoInterp );
	//	        }

	//	        // assign the primary color (PER_VERTEX required for later optimization)
	//	        osg::Vec4Array* colors = new osg::Vec4Array;
	//	        colors->assign( osgGeom->getVertexArray()->getNumElements(), primaryColor );
	//	        osgGeom->setColorArray( colors );
	//	        osgGeom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );

	//	        geode->addDrawable( osgGeom );

	//	        // record the geometry's primitive set(s) in the index:
	//	        if ( context.featureIndex() )
	//	            context.featureIndex()->tagDrawable( osgGeom, input );
	//	    
	//	        // install clamping attributes if necessary
	//	        if (_style.has<AltitudeSymbol>() &&
	//	            _style.get<AltitudeSymbol>()->technique() == AltitudeSymbol::TECHNIQUE_GPU)
	//	        {
	//	            Clamping::applyDefaultClampingAttrs( osgGeom, input->getDouble("__oe_verticalOffset", 0.0) );
	//	        }
	//	    }
	//	}
	//}
 //   return geode;

}

ZoomUpdateCallback* BuildGeometryFilter::_zoomUpdateCallback = NULL;
float ZoomUpdateCallback::_scale=1.0;
float  BuildGeometryFilter::_hscale = 1.0;

SineAnimation* BuildGeometryFilter::_sineAnimation = NULL;
float SineAnimation::_stepsize=1.0;
float  BuildGeometryFilter::_hstepsize = 1.0;

osg::Geode*
BuildGeometryFilter::processPoints(FeatureList& features, FilterContext& context)
{
    osg::Geode* geode = new osg::Geode();

    bool makeECEF = false;
    const SpatialReference* featureSRS = 0L;
    const SpatialReference* mapSRS = 0L;

    // set up referencing information:
    if ( context.isGeoreferenced() )
    {
        makeECEF   = context.getSession()->getMapInfo().isGeocentric();
        featureSRS = context.extent()->getSRS();
        mapSRS     = context.getSession()->getMapInfo().getProfile()->getSRS();
    }
#pragma region 正常
	if(!_style.get<PointSymbol>()->isMutliple().get()&&!_style.get<PointSymbol>()->isEfficient().get())
	{
		for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
		{
			Feature* input = f->get();

			GeometryIterator parts( input->getGeometry(), true );
			while( parts.hasMore() )
			{
				Geometry* part = parts.next();

				// extract the required point symbol; bail out if not found.
				const PointSymbol* point =
					input->style().isSet() && input->style()->has<PointSymbol>() ? input->style()->get<PointSymbol>() :
					_style.get<PointSymbol>();

				if ( !point )
					continue;

				// resolve the color:
				osg::Vec4f primaryColor = point->fill()->color();
            
				osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
				osgGeom->setUseVertexBufferObjects( true );
				osgGeom->setUseDisplayList( false );

				// embed the feature name if requested. Warning: blocks geometry merge optimization!
				if ( _featureNameExpr.isSet() )
				{
					const std::string& name = input->eval( _featureNameExpr.mutable_value(), &context );
					osgGeom->setName( name );
				}

				// build the geometry:
				osg::Vec3Array* allPoints = new osg::Vec3Array();

				transformAndLocalize( part->asVector(), featureSRS, allPoints, mapSRS, _world2local, makeECEF );

				osgGeom->addPrimitiveSet( new osg::DrawArrays(GL_POINTS, 0, allPoints->getNumElements()) );
				osgGeom->setVertexArray( allPoints );

				if ( input->style().isSet() )
				{
					//TODO: re-evaluate this. does it hinder geometry merging?
					applyPointSymbology( osgGeom->getOrCreateStateSet(), point );
				}

				// assign the primary color (PER_VERTEX required for later optimization)
				osg::Vec4Array* colors = new osg::Vec4Array;
				colors->assign( osgGeom->getVertexArray()->getNumElements(), primaryColor );
				osgGeom->setColorArray( colors );
				osgGeom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );

				geode->addDrawable( osgGeom );

				// record the geometry's primitive set(s) in the index:
				if ( context.featureIndex() )
					context.featureIndex()->tagDrawable( osgGeom, input );
        
				// install clamping attributes if necessary
				if (_style.has<AltitudeSymbol>() &&
					_style.get<AltitudeSymbol>()->technique() == AltitudeSymbol::TECHNIQUE_GPU)
				{            
					Clamping::applyDefaultClampingAttrs( osgGeom, input->getDouble("__oe_verticalOffset", 0.0) );
				}
			}
		}
	}
#pragma endregion
#pragma region 高效化渲染
	else if(!_style.get<PointSymbol>()->isMutliple().get()&&_style.get<PointSymbol>()->isEfficient().get())
	{
		osg::Vec3Array* allPoints = new osg::Vec3Array();
		std::vector<osg::Vec3d> pointsCon;
		const PointSymbol* point = _style.get<PointSymbol>();
		osg::Vec4f primaryColor = point->fill()->color();
		osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
		osgGeom->setUseVertexBufferObjects( true );
		osgGeom->setUseDisplayList( false );


		for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
		{
			Feature* input = f->get();
			GeometryIterator parts( input->getGeometry(), true );
			while( parts.hasMore() )
			{
				Geometry* part = parts.next();
				std::copy( part->asVector().begin(), part->asVector().end(),std::back_inserter(pointsCon));
			}
		}
		transformAndLocalize( pointsCon, featureSRS, allPoints, mapSRS, _world2local, makeECEF );
		osgGeom->addPrimitiveSet( new osg::DrawArrays(GL_POINTS, 0, allPoints->getNumElements()) );
		osgGeom->setVertexArray( allPoints );

		osg::Vec4Array* colors = new osg::Vec4Array;
		colors->assign( osgGeom->getVertexArray()->getNumElements(), primaryColor );
		osgGeom->setColorArray( colors );
		osgGeom->setColorBinding( osg::Geometry::BIND_PER_VERTEX );
		geode->addDrawable( osgGeom );
	}
	//流场动态可视化
	else if (_style.get<PointSymbol>()->isFlow().get())
	{
		osg::Vec3Array* allPoints = new osg::Vec3Array();
		std::vector<osg::Vec3d> pointsCon;
		const PointSymbol* point = _style.get<PointSymbol>();
		
		
		for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
		{
			Feature* input = f->get();
			GeometryIterator parts( input->getGeometry(), true );
			while( parts.hasMore() )
			{
				Geometry* part = parts.next();	
				std::copy( part->asVector().begin(), part->asVector().end(),std::back_inserter(pointsCon));
			}
		}

		transformAndLocalize( pointsCon, featureSRS, allPoints, mapSRS, _world2local, makeECEF );
	
		osg::Vec4Array* colorArray = new osg::Vec4Array;
		osg::FloatArray* speedArray = new osg::FloatArray;
		int speedStep = 0;
		int count = allPoints->getTotalDataSize()/allPoints->getDataSize()/4;         
		for (int i=0;i<count;i++)
		{
			osg::Vec3f point1 = allPoints->at(i);
			point1 = point1 + osg::Vec3f(0.0,0.0,1.5);
			allPoints->at(i) = point1; 

			//speed
			int speedStart =  point->rotationString().get().find('*',speedStep);
			int speedEnd =  point->rotationString().get().find('*',speedStart+1);
			speedStep = speedEnd;

			float speed;
			std::string tempSpeedCell = point->rotationString().get().substr(speedStart+1,speedEnd-speedStart-1);
			std::istringstream iss1(tempSpeedCell);
			iss1>>speed;
			speedArray->push_back(speed);
		}

		int sizeStep = 0;
		double geomcount;
		osg::FloatArray* sizeArray = new osg::FloatArray;
		
		if (point->isFeatureNode().get())
		{
			geomcount = point->featureCount().get();
		}
		else
		{
			geomcount = features.size();
		}
		for(int k = 0;k<geomcount;k++)
		{
			int sizeStart =  point->sizeString().get().find('*',sizeStep);
			int sizeEnd =  point->sizeString().get().find('*',sizeStart+1);
			sizeStep = sizeEnd;

			float size;
			std::string tempSizeCell = point->sizeString().get().substr(sizeStart+1,sizeEnd-sizeStart-1);
			std::istringstream iss1(tempSizeCell);
			iss1>>size;
			sizeArray->push_back(size);
		}

		osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
		osgGeom->setUseVertexBufferObjects( true );
		osgGeom->setUseDisplayList( false );
		
		int geomindex = 0;
		int localcount = 0;
		int k = 0;
		int m = 0;
		for (int i=0;i<sizeArray->size();i++)
		{
			double realcount = 0;
			double nowcount = 0;
			for (int j=0;j<sizeArray->at(i);j++)
			{
				realcount += 1.0/speedArray->at(k);
				k++;
			}

			for (int j=0;j<sizeArray->at(i);j++)
			{
				nowcount += 1.0/speedArray->at(m);
				
				//速度据观察都在0-1
				if (speedArray->at(m)>1.0)
				{
					//osg::Vec4 trueColor = osg::Vec4(1.0/sizeArray->at(i),1.0,0.0,1.0/(1+j));
					osg::Vec4 trueColor = osg::Vec4(1.0/realcount,1.0,0.0,1.0/(1+nowcount));
					colorArray->push_back(trueColor);
				}
				else
				{
					//osg::Vec4 trueColor = osg::Vec4(1.0/sizeArray->at(i),speedArray->at(geomindex+j),0.0,1.0/(1+j));
					osg::Vec4 trueColor = osg::Vec4(1.0/realcount,speedArray->at(m),0.0,1.0/(1+nowcount));
					colorArray->push_back(trueColor);
				}
				m++;
			}
			osgGeom->addPrimitiveSet( new osg::DrawArrays(GL_LINE_STRIP, geomindex, sizeArray->at(i)) );
			geomindex += sizeArray->at(i);
			//geomindex += realcount;
		}


		osgGeom->setVertexArray( allPoints );
		osg::StateSet *ss = osgGeom->getOrCreateStateSet();

		//统一一致变量回调
		osg::Uniform* u_anim1( new osg::Uniform( "u_anim1", 0.0f ) );
		//u_anim1->setUpdateCallback( new SineAnimation( 0.0f ) );

		_sineAnimation = new SineAnimation();
		u_anim1->setUpdateCallback(_sineAnimation);

		ss->addUniform( u_anim1 );
		ss->setAttributeAndModes(new osg::LineWidth(2.0f), 1);
		
		osg::Program *program = new osg::Program;
		char* fragShader;
		char* vertShader;
		char* geoShader;
		
		program->addBindAttribLocation("a_color",1);
		osgGeom->setVertexAttribArray(1,colorArray);
		osgGeom->setVertexAttribBinding(1,osg::Geometry::BIND_PER_VERTEX);

		fragShader = 
			"#version 330 compatibility\n"
			"varying vec4 v_color_out;\n"
			"void main(void)\n"
			"{\n"
			"	gl_FragColor = gl_Color;\n"
			"}\n"            
			;
		vertShader = 
			"#version 330 compatibility\n"
			"varying vec4 v_color;\n"
			"attribute vec4 a_color;\n"
			"void main(void)\n"
			"{\n" 
			"   gl_Position = gl_Vertex;\n"
			"	gl_FrontColor = a_color;\n"
			"	gl_TexCoord[0] = a_color;\n"
			"}\n"
			;
		geoShader = 
			"#version 330 compatibility\n"
			"#extension GL_EXT_geometry_shader4 : enable\n"
			"uniform float u_anim1;\n"
			"layout(lines) in;\n"
			"layout(line_strip,max_vertices=1024) out;\n"
			"varying in vec4 v_color[];\n"
			"varying out vec4 v_color_out;\n"
			"void main(void)\n"
			"{\n"
			"    vec4 v = gl_PositionIn[0];\n"
			"	int j = int(u_anim1-0.2)%(int(1.0/gl_TexCoordIn[0][0].x)+20);\n"
			"	int index = int(1.0/gl_TexCoordIn[0][0].w-1);\n"
			"	for (int i=0;i<gl_VerticesIn;i++)\n"
				"{\n"
				"	gl_Position = gl_ModelViewProjectionMatrix * gl_PositionIn[i];\n"
				"if (index<j)\n"
				"{\n"
				/*"if(1.0/gl_TexCoordIn[0][0].x>500)\n"
				"{\n"
				"		gl_FrontColor = vec4(1.0,0.0,0.0,((index+1.0)*sqrt(abs(gl_TexCoordIn[0][0].y)))/j);\n"
				"}\n"
				"else\n"
				"{\n"*/
				"		gl_FrontColor = vec4(1.0,1.0,1.0,((index+1.0)*sqrt(abs(gl_TexCoordIn[0][0].y)))/j);\n"
				/*"}\n"*/
				"}\n"
				"else\n" 
				"{\n"
				"		gl_FrontColor = vec4(1.0,1.0,1.0,0.0);\n"
				"}\n"
				"    EmitVertex();\n"
				"}\n"
				"   EndPrimitive();\n" 
			"}\n"
			;
		ss->setMode(GL_LINE_SMOOTH,osg::StateAttribute::ON);
		program->addShader(new osg::Shader(osg::Shader::FRAGMENT,fragShader));
		program->addShader(new osg::Shader(osg::Shader::VERTEX,vertShader));
		program->addShader(new osg::Shader(osg::Shader::GEOMETRY,geoShader));
		ss->setAttributeAndModes(program,osg::StateAttribute::ON);
		
		geode->addDrawable( osgGeom );
	}
#pragma endregion
#pragma region 符号化
	else
	{
		osg::Vec3Array* allPoints = new osg::Vec3Array();
		std::vector<osg::Vec3d> pointsCon;
		const PointSymbol* point = _style.get<PointSymbol>();
		
		for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
		{
			Feature* input = f->get();
			GeometryIterator parts( input->getGeometry(), true );
			while( parts.hasMore() )
			{
				Geometry* part = parts.next();
				std::copy( part->asVector().begin(), part->asVector().end(),std::back_inserter(pointsCon));
			}
		}
		pointsCon.push_back(osg::Vec3d(120,30,-6371000));
		transformAndLocalize( pointsCon, featureSRS, allPoints, mapSRS, _world2local, makeECEF );
	
		osg::Vec3f earthHeart = allPoints->at(allPoints->size()-1);
		allPoints->pop_back();

		int count = allPoints->getTotalDataSize()/allPoints->getDataSize()/4;      //add by jxf 为解决被遮挡问题    
		for (int i=0;i<count;i++)
		{
			osg::Vec3f point = allPoints->at(i);
			point = point + osg::Vec3f(0.0,0.0,1.0);
			allPoints->at(i) = point; 
		}

		osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
		osgGeom->setUseVertexBufferObjects( true );
		osgGeom->setUseDisplayList( false );
		osgGeom->addPrimitiveSet( new osg::DrawArrays(GL_POINTS, 0, allPoints->getNumElements()) );
		osgGeom->setVertexArray( allPoints );

		
		osg::StateSet *ss = osgGeom->getOrCreateStateSet();
		if (/*point->scale().get()*/false)
		{
			/*GLuint ShaderProgram = glCreateProgram();
			glUseProgram(ShaderProgram);
			GLint VertexColorLocation = glGetUniformLocation(ShaderProgram,"scale");
			float scale = point->scale().get();
			glUniform1f(VertexColorLocation,scale);*/
		/*	osg::Uniform* uni = ss->getUniform("scale");
			osg::Uniform* uni1 = ss->getOrCreateUniform("scale");
			int i = ss->getUniformList().size();
			
			float scale = point->scale().get();
			uni->set(scale);*/
			//int i = ss->getUniformList().size();
			//_zoomUpdateCallback=(ZoomUpdateCallback*)(ss->getUniform("scale")->getUpdateCallback());
			//_zoomUpdateCallback->_scale= point->scale().get();
		}
		else
		{
			osg::Program *program = new osg::Program;
			char* fragShader;
			char* vertShader;
			char* geoShader;
		
			//ss->setAttributeAndModes(new osg::Point(10), 1);
			//osg::Vec4f mainColor = point->fill()->color();
			int colorStep = 0;int rotationStep = 0;int sizeStep = 0;int markerTypeStep = 0;
			osg::Vec4Array* colorArray = new osg::Vec4Array;
			osg::FloatArray* sizeArray = new osg::FloatArray;
			osg::FloatArray* rotationArray = new osg::FloatArray;
			osg::FloatArray* markerTypeArray = new osg::FloatArray;
			double size;
			if (point->isFeatureNode().get())
			{
				size = point->featureCount().get();
			}
			else
			{
				size = features.size();
			}
			for(int k = 0;k<size;k++)
			{
				int colorStart = point->colorString().get().find('*',colorStep);
				int colorEnd = point->colorString().get().find('*',colorStart+1);
				colorStep = colorEnd;
				int sizeStart =  point->sizeString().get().find('*',sizeStep);
				int sizeEnd =  point->sizeString().get().find('*',sizeStart+1);
				sizeStep = sizeEnd;
				int rotationStart =  point->rotationString().get().find('*',rotationStep);
				int rotationEnd =  point->rotationString().get().find('*',rotationStart+1);
				rotationStep = rotationEnd;
				int markerTypeStart =  point->markerTypeString().get().find('*',markerTypeStep);
				int markerTypeEnd =  point->markerTypeString().get().find('*',markerTypeStep+1);
				markerTypeStep = markerTypeEnd;

				float size,rotation,markerTypeStep;
				std::string tempSizeCell = point->sizeString().get().substr(sizeStart+1,sizeEnd-sizeStart-1);
				std::istringstream iss1(tempSizeCell);
				iss1>>size;
				//std::cout<<size<<std::endl;
				sizeArray->push_back(size);
				std::string tempRotationCell = point->rotationString().get().substr(rotationStart+1,rotationEnd-rotationStart-1);
				std::istringstream iss2(tempRotationCell);
				iss2>>rotation;
				//std::cout<<rotation<<std::endl;
				rotationArray->push_back(rotation);
				std::string tempMarkerTypeCell = point->markerTypeString().get().substr(markerTypeStart+1,markerTypeEnd-markerTypeStart-1);
				std::istringstream iss3(tempMarkerTypeCell);
				iss3>>markerTypeStep;
				//std::cout<<markerTypeStep<<std::endl;
				markerTypeArray->push_back(markerTypeStep);

				float trueColorFloat[3];
				std::string tempVecString = point->colorString().get().substr(colorStart+1,colorEnd-colorStart-1);
				{
					int tempStep = 0;

					for(int i=0;i<3;i++)
					{
						int tempStart = tempVecString.find(',',tempStep);
						int tempEnd = tempVecString.find(',',tempStart+1);
						tempStep = tempEnd;
						std::string tempColorCell = tempVecString.substr(tempStart+1,tempEnd-tempStart-1);
						std::istringstream iss(tempColorCell);
						iss>>trueColorFloat[i];
					}
				}

				osg::Vec4 trueColor = osg::Vec4(trueColorFloat[0],trueColorFloat[1],trueColorFloat[2],1.0);
				colorArray->push_back(trueColor);
			}

			program->addBindAttribLocation("a_color",1);
			osgGeom->setVertexAttribArray(1,colorArray);
			osgGeom->setVertexAttribBinding(1,osg::Geometry::BIND_PER_VERTEX);

			program->addBindAttribLocation("a_size",2);
			osgGeom->setVertexAttribArray(2,sizeArray);
			osgGeom->setVertexAttribBinding(2,osg::Geometry::BIND_PER_VERTEX);

			program->addBindAttribLocation("a_rotation",3);
			osgGeom->setVertexAttribArray(3,rotationArray);
			osgGeom->setVertexAttribBinding(3,osg::Geometry::BIND_PER_VERTEX);

			program->addBindAttribLocation("a_markerType",4);
			osgGeom->setVertexAttribArray(4,markerTypeArray);
			osgGeom->setVertexAttribBinding(4,osg::Geometry::BIND_PER_VERTEX);


			osg::Uniform* uearthHeart = new osg::Uniform("earthHeart",earthHeart);
			ss->addUniform(uearthHeart);
			float scale = 1.0;
			if (point->scale().get())
			{
				scale = point->scale().get();
			}
			osg::Uniform* uscale = new osg::Uniform("scale",scale);
			_zoomUpdateCallback = new ZoomUpdateCallback();
			uscale->setUpdateCallback(_zoomUpdateCallback);
			ss->addUniform(uscale);
		
			fragShader = 
				"#version 330 compatibility\n"
				//"uniform int symbolType;\n"
				"void main(void){\n"
				"	float markerType = gl_TexCoord[0].z;\n"
				//"	gl_FragDepth = clamp(gl_TexCoord[0].w,0.0,1.0);\n"
				//"	if(gl_TexCoord[0].w>gl_FragDepth)\n"
				//"		discard;\n"
				"	if(markerType==0.0){\n"
				"		vec3 N;\n"
				"		N.xy = gl_TexCoord[0].xy*vec2(2.0, -2.0) + vec2(-1.0, 1.0);\n"
				"		float r2 = dot(N.xy, N.xy);\n"
				"		if (r2 > 1.0) discard;\n"
				"	}\n"
				"	if(markerType>=1.5&&markerType<2.5){\n"
				"		if ((0.5-gl_TexCoord[0].x)*2-gl_TexCoord[0].y>=0 || (gl_TexCoord[0].x-0.5)*2-gl_TexCoord[0].y>=0)\n"
				"			discard;\n"
				"	}\n"
				"	if(markerType>=2.5&&markerType<3.5){\n"
				"		float Angle = radians(72.0);\n"// / 180 * 3.1415926535897932;
				"		float Ax = 0.5;\n"
				"		float Ay = 1.0;\n"
				"		float Bx =  (Ax - 0.5) * cos(Angle) + (Ay - 0.5) * sin(Angle) + 0.5;\n"
				"		float By = -(Ax - 0.5) * sin(Angle) + (Ay - 0.5) * cos(Angle) + 0.5;\n"
				"		float Cx =  (Bx - 0.5) * cos(Angle) + (By - 0.5) * sin(Angle) + 0.5;\n"
				"		float Cy = -(Bx - 0.5) * sin(Angle) + (By - 0.5) * cos(Angle) + 0.5;\n"
				"		float Dx =  (Cx - 0.5) * cos(Angle) + (Cy - 0.5) * sin(Angle) + 0.5;\n"
				"		float Dy = -(Cx - 0.5) * sin(Angle) + (Cy - 0.5) * cos(Angle) + 0.5;\n"
				"		float Ex =  (Dx - 0.5) * cos(Angle) + (Dy - 0.5) * sin(Angle) + 0.5;\n"
				"		float Ey = -(Dx - 0.5) * sin(Angle) + (Dy - 0.5) * cos(Angle) + 0.5;\n"

				"		float kDA = (Dy - Ay) * 1.0 / (Dx - Ax);\n"
				"		float bDA = Dy - kDA * Dx;\n"

				"		float kAC = (Ay - Cy) * 1.0 / (Ax - Cx);\n"
				"		float bAC = Ay - kAC * Ax;\n"

				"		float kEB = (Ey - By) * 1.0 / (Ex - Bx);\n"
				"		float bEB = Ey - kEB * Ex;\n"

				"		float kBD = (By - Dy) * 1.0 / (Bx - Dx);\n"
				"		float bBD = By - kBD * Bx;\n"

				"		float kEC = (Ey - Cy) * 1.0 / (Ex - Cx);\n"
				"		float bEC = Ey - kEC * Ex;\n"

				"		int flag = 0;\n"
				"		float x = gl_TexCoord[0].x;\n"
				"		float y = gl_TexCoord[0].y;\n"
				"		if (y <= kAC * x + bAC && y <= kDA * x + bDA && (y >= kBD * x + bBD || y >= kEC * x + bEC))\n"
				"		{\n"
				"			flag = 1;\n"
				"		}\n"
				"		else if (y <= kEB * x + bEB)\n"
				"		{\n"
				"			if ((y >= kEC * x + bEC && y >= kDA * x + bDA) || (y >= kAC * x + bAC && y >= kBD * x + bBD))\n"
				"			{\n"
				"				flag = 1;\n"
				"			}\n"
				"			else\n"
				"			{\n"
				"				flag = 0;\n"
				"			}\n"
				"		}\n"
				"		else\n"
				"		{\n"
				"			flag = 0;\n"
				"		}\n"
				"		if(flag<1)\n"
				"			discard;\n"
				"	}\n"

				"	if(markerType>=3.5&&markerType<4.5){\n"
				"		if ((0.5-gl_TexCoord[0].x)*2-gl_TexCoord[0].y>=0 || (gl_TexCoord[0].x-0.5)*2-gl_TexCoord[0].y>=0)\n"
				"			discard;\n"
				"		if (gl_TexCoord[0].x<=0.5&&gl_TexCoord[0].y+2.0/3*gl_TexCoord[0].x-1>=0)\n"
				"			discard;\n"
				"		if (gl_TexCoord[0].x>=0.5&&gl_TexCoord[0].y-2.0/3*gl_TexCoord[0].x-1.0/3>=0)\n"
				"			discard;\n"
				"	}\n"
				"	gl_FragColor = gl_Color;\n"
				"}\n"
				;
			vertShader = 
				"#version 330 compatibility\n"
				"attribute vec4 a_color;\n"
				"attribute float a_size;\n"
				"attribute float a_rotation;\n"
				"attribute float a_markerType;\n"
				//"varying float v_size;\n"
				//"varying float v_rotation;\n"
				"void main(){\n"
				//"	v_size = a_size;\n"
				//"	v_rotation = a_rotation;\n"
				"	gl_TexCoord[0] = vec4(a_size,a_rotation,a_markerType, 1.0);\n"
				"	gl_FrontColor = a_color;\n"
				"	gl_Position = gl_Vertex;\n"
				"}\n"
				;
			geoShader = 
				"#version 330 compatibility\n"
				"#extension GL_EXT_geometry_shader4:enable\n"
				"layout(points) in;\n"
				"layout(triangle_strip,max_vertices=4) out;\n"
				"uniform float scale;"
				"uniform vec3 earthHeart;"
				"void main(){\n"
				//"	testVec = gl_ModelViewMatrix * gl_ProjectionMatrix * testVec;\n"
				"	float r = sqrt(earthHeart.x*earthHeart.x+earthHeart.y*earthHeart.y+earthHeart.z*earthHeart.z);\n"
				"	float a = gl_PositionIn[0].x-earthHeart.x;\n"
				"	float b = gl_PositionIn[0].y-earthHeart.y;\n"
				"	float c = gl_PositionIn[0].z-earthHeart.z;\n"
				"	float w = gl_PositionIn[0].w;\n"

				"	float s = gl_TexCoordIn[0][0].x/30000*scale;\n"
				"	float ra = radians(gl_TexCoordIn[0][0].y-45);\n"
				"	float rb = radians(gl_TexCoordIn[0][0].y+45);\n"
				"	float rc = radians(gl_TexCoordIn[0][0].y+135);\n"
				"	float rd = radians(gl_TexCoordIn[0][0].y+225);\n"
				"	float markerType = gl_TexCoordIn[0][0].z;\n"

				"	if((r*r)/b-b>0)\n"
				"	{\n"
				"		gl_TexCoord[0] = vec4(1, 1, markerType, 1);\n "
				"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rc)+r*c/b*sin(rc)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rc),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rc)-r*a/b*sin(rc)),w);\n"
				"		gl_FrontColor = gl_FrontColorIn[0];\n"  
				"		EmitVertex();\n" 

				"		gl_TexCoord[0] = vec4(1, 0, markerType, 1);\n "
				"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rb)+r*c/b*sin(rb)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rb),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rb)-r*a/b*sin(rb)),w);\n;"
				"		gl_FrontColor = gl_FrontColorIn[0];\n"  
				"		EmitVertex();\n"

				"		gl_TexCoord[0] = vec4(0, 1, markerType, 1);\n "
				"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rd)+r*c/b*sin(rd)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rd),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rd)-r*a/b*sin(rd)),w);\n"
				"		gl_FrontColor = gl_FrontColorIn[0];\n"  
				"		EmitVertex();\n"

				"		gl_TexCoord[0] = vec4(0, 0, markerType, 1);\n "
				"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(ra)+r*c/b*sin(ra)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(ra),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(ra)-r*a/b*sin(ra)),w);\n"
				"		gl_FrontColor = gl_FrontColorIn[0];\n" 
				"		EmitVertex();\n"
				"	}\n"

				"	else\n"
				"	{\n"
				"		gl_TexCoord[0] = vec4(1, 0, markerType, 1);\n "
				"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rc)+r*c/b*sin(rc)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rc),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rc)-r*a/b*sin(rc)),w);\n"
				"		gl_FrontColor = gl_FrontColorIn[0];\n"  
				"		EmitVertex();\n"

				"		gl_TexCoord[0] = vec4(1, 1, markerType, 1);\n "
				"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rb)+r*c/b*sin(rb)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rb),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rb)-r*a/b*sin(rb)),w);\n;"
				"		gl_FrontColor = gl_FrontColorIn[0];\n"  
				"		EmitVertex();\n"

				"		gl_TexCoord[0] = vec4(0, 0, markerType, 1);\n "
				"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rd)+r*c/b*sin(rd)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rd),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rd)-r*a/b*sin(rd)),w);\n"
				"		gl_FrontColor = gl_FrontColorIn[0];\n"  
				"		EmitVertex();\n"

				"		gl_TexCoord[0] = vec4(0, 1, markerType, 1);\n "
				"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(ra)+r*c/b*sin(ra)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(ra),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(ra)-r*a/b*sin(ra)),w);\n"
				"		gl_FrontColor = gl_FrontColorIn[0];\n"  
				"		EmitVertex();\n"
				"	}\n"
				"	EndPrimitive();\n"
				"}\n"
				;


			program->addShader(new osg::Shader(osg::Shader::FRAGMENT,fragShader));
			program->addShader(new osg::Shader(osg::Shader::VERTEX,vertShader));
			program->addShader(new osg::Shader(osg::Shader::GEOMETRY,geoShader));
			ss->setAttributeAndModes(program,osg::StateAttribute::ON);
		}

		geode->addDrawable( osgGeom );
	}
//else
//	{
//		int i =0;
//		for( FeatureList::iterator f = features.begin(); f != features.end(); ++f )
//		{
//			osg::Vec3Array* allPoints = new osg::Vec3Array();
//			std::vector<osg::Vec3d> pointsCon;
//			const PointSymbol* point = _style.get<PointSymbol>();
//			Feature* input = f->get();
//			GeometryIterator parts( input->getGeometry(), true );
//			while( parts.hasMore() )
//			{
//				Geometry* part = parts.next();
//				std::copy( part->asVector().begin(), part->asVector().end(),std::back_inserter(pointsCon));
//				pointsCon.push_back(osg::Vec3d(120,30,-6371000));
//				transformAndLocalize( pointsCon, featureSRS, allPoints, mapSRS, _world2local, makeECEF );
//				
//				osg::Vec3f earthHeart = allPoints->at(allPoints->size()-1);
//				allPoints->pop_back();
//				int count = allPoints->getTotalDataSize()/allPoints->getDataSize()/4;      //add by jxf 为解决被遮挡问题    
//				for (int i=0;i<count;i++)
//				{
//					osg::Vec3f point = allPoints->at(i);
//					point = point + osg::Vec3f(0.0,0.0,2.0);
//					allPoints->at(i) = point; 
//				}
//				osg::ref_ptr<osg::Geometry> osgGeom = new osg::Geometry();
//				osgGeom->setUseVertexBufferObjects( true );
//				osgGeom->setUseDisplayList( false );
//
//				if ( _featureNameExpr.isSet() )
//				{
//					const std::string& name = input->eval( _featureNameExpr.mutable_value(), &context );
//					osgGeom->setName( name );
//				}
//			
//				osgGeom->addPrimitiveSet( new osg::DrawArrays(GL_POINTS, 0, allPoints->getNumElements()) );
//				osgGeom->setVertexArray( allPoints );
//
//				if ( input->style().isSet() )
//				{
//					//TODO: re-evaluate this. does it hinder geometry merging?
//					applyPointSymbology( osgGeom->getOrCreateStateSet(), point );
//				}
//
//				osg::StateSet *ss = osgGeom->getOrCreateStateSet();
//				//ss->setAttributeAndModes(new osg::Point(10), 1);
//				osg::Program *program = new osg::Program;
//		
//				//osg::Vec4f mainColor = point->fill()->color();
//				int colorStep = 0;int rotationStep = 0;int sizeStep = 0;int markerTypeStep = 0;
//				osg::Vec4Array* colorArray = new osg::Vec4Array;
//				osg::FloatArray* sizeArray = new osg::FloatArray;
//				osg::FloatArray* rotationArray = new osg::FloatArray;
//				osg::FloatArray* markerTypeArray = new osg::FloatArray;
//				for(int k = 0;k<features.size();k++)
//				{
//					int colorStart = point->colorString().get().find('*',colorStep);
//					int colorEnd = point->colorString().get().find('*',colorStart+1);
//					colorStep = colorEnd;
//					int sizeStart =  point->sizeString().get().find('*',sizeStep);
//					int sizeEnd =  point->sizeString().get().find('*',sizeStart+1);
//					sizeStep = sizeEnd;
//					int rotationStart =  point->rotationString().get().find('*',rotationStep);
//					int rotationEnd =  point->rotationString().get().find('*',rotationStart+1);
//					rotationStep = rotationEnd;
//					int markerTypeStart =  point->markerTypeString().get().find('*',markerTypeStep);
//					int markerTypeEnd =  point->markerTypeString().get().find('*',markerTypeStep+1);
//					markerTypeStep = markerTypeEnd;
//
//					if(k == i)
//					{
//						float size,rotation,markerTypeStep;
//						std::string tempSizeCell = point->sizeString().get().substr(sizeStart+1,sizeEnd-sizeStart-1);
//						std::istringstream iss1(tempSizeCell);
//						iss1>>size;
//						//std::cout<<size<<std::endl;
//						sizeArray->push_back(size);
//						std::string tempRotationCell = point->rotationString().get().substr(rotationStart+1,rotationEnd-rotationStart-1);
//						std::istringstream iss2(tempRotationCell);
//						iss2>>rotation;
//						//std::cout<<rotation<<std::endl;
//						rotationArray->push_back(rotation);
//						std::string tempMarkerTypeCell = point->markerTypeString().get().substr(markerTypeStart+1,markerTypeEnd-markerTypeStart-1);
//						std::istringstream iss3(tempMarkerTypeCell);
//						iss3>>markerTypeStep;
//						//std::cout<<markerTypeStep<<std::endl;
//						markerTypeArray->push_back(markerTypeStep);
//
//						float trueColorFloat[3];
//						std::string tempVecString = point->colorString().get().substr(colorStart+1,colorEnd-colorStart-1);
//						{
//							int tempStep = 0;
//
//							for(int i=0;i<3;i++)
//							{
//								int tempStart = tempVecString.find(',',tempStep);
//								int tempEnd = tempVecString.find(',',tempStart+1);
//								tempStep = tempEnd;
//								std::string tempColorCell = tempVecString.substr(tempStart+1,tempEnd-tempStart-1);
//								std::istringstream iss(tempColorCell);
//								iss>>trueColorFloat[i];
//							}
//						}
//
//						osg::Vec4 trueColor = osg::Vec4(trueColorFloat[0],trueColorFloat[1],trueColorFloat[2],1.0);
//						colorArray->push_back(trueColor);
//						break;
//					}
//				}
//				i++;
//				program->addBindAttribLocation("a_color",1);
//				osgGeom->setVertexAttribArray(1,colorArray);
//				osgGeom->setVertexAttribBinding(1,osg::Geometry::BIND_PER_VERTEX);
//
//				program->addBindAttribLocation("a_size",2);
//				osgGeom->setVertexAttribArray(2,sizeArray);
//				osgGeom->setVertexAttribBinding(2,osg::Geometry::BIND_PER_VERTEX);
//
//				program->addBindAttribLocation("a_rotation",3);
//				osgGeom->setVertexAttribArray(3,rotationArray);
//				osgGeom->setVertexAttribBinding(3,osg::Geometry::BIND_PER_VERTEX);
//
//				program->addBindAttribLocation("a_markerType",4);
//				osgGeom->setVertexAttribArray(4,markerTypeArray);
//				osgGeom->setVertexAttribBinding(4,osg::Geometry::BIND_PER_VERTEX);
//
//
//				osg::Uniform* uearthHeart = new osg::Uniform("earthHeart",earthHeart);
//				ss->addUniform(uearthHeart);
//				char* fragShader;
//				char* vertShader;
//				char* geoShader;
//					fragShader = 
//						"#version 330 compatibility\n"
//						//"uniform int symbolType;\n"
//						"flat in vec4 oe_pick_encoded_objectid; \n"
//						"flat in int  oe_pick_color_contains_objectid; \n"
//
//						"out vec4 fragColor; \n"
//						"void main(void){\n"
//
//						"    if ( oe_pick_color_contains_objectid == 1 ) \n"
//						"        fragColor = gl_Color; \n"
//						"    else \n"
//						"        fragColor = oe_pick_encoded_objectid; \n"
//
//						"	float markerType = gl_TexCoord[0].z;\n"
//						//"	gl_FragDepth = clamp(gl_TexCoord[0].w,0.0,1.0);\n"
//						//"	if(gl_TexCoord[0].w>gl_FragDepth)\n"
//						//"		discard;\n"
//						"	if(markerType==0.0){\n"
//						"		vec3 N;\n"
//						"		N.xy = gl_TexCoord[0].xy*vec2(2.0, -2.0) + vec2(-1.0, 1.0);\n"
//						"		float r2 = dot(N.xy, N.xy);\n"
//						"		if (r2 > 1.0) discard;\n"
//						"	}\n"
//						"	if(markerType>=1.5&&markerType<2.5){\n"
//						"		if ((0.5-gl_TexCoord[0].x)*2-gl_TexCoord[0].y>=0 || (gl_TexCoord[0].x-0.5)*2-gl_TexCoord[0].y>=0)\n"
//						"			discard;\n"
//						"	}\n"
//						"	if(markerType>=2.5&&markerType<3.5){\n"
//						"		float Angle = radians(72.0);\n"// / 180 * 3.1415926535897932;
//						"		float Ax = 0.5;\n"
//						"		float Ay = 1.0;\n"
//						"		float Bx =  (Ax - 0.5) * cos(Angle) + (Ay - 0.5) * sin(Angle) + 0.5;\n"
//						"		float By = -(Ax - 0.5) * sin(Angle) + (Ay - 0.5) * cos(Angle) + 0.5;\n"
//						"		float Cx =  (Bx - 0.5) * cos(Angle) + (By - 0.5) * sin(Angle) + 0.5;\n"
//						"		float Cy = -(Bx - 0.5) * sin(Angle) + (By - 0.5) * cos(Angle) + 0.5;\n"
//						"		float Dx =  (Cx - 0.5) * cos(Angle) + (Cy - 0.5) * sin(Angle) + 0.5;\n"
//						"		float Dy = -(Cx - 0.5) * sin(Angle) + (Cy - 0.5) * cos(Angle) + 0.5;\n"
//						"		float Ex =  (Dx - 0.5) * cos(Angle) + (Dy - 0.5) * sin(Angle) + 0.5;\n"
//						"		float Ey = -(Dx - 0.5) * sin(Angle) + (Dy - 0.5) * cos(Angle) + 0.5;\n"
//
//						"		float kDA = (Dy - Ay) * 1.0 / (Dx - Ax);\n"
//						"		float bDA = Dy - kDA * Dx;\n"
//
//						"		float kAC = (Ay - Cy) * 1.0 / (Ax - Cx);\n"
//						"		float bAC = Ay - kAC * Ax;\n"
//
//						"		float kEB = (Ey - By) * 1.0 / (Ex - Bx);\n"
//						"		float bEB = Ey - kEB * Ex;\n"
//
//						"		float kBD = (By - Dy) * 1.0 / (Bx - Dx);\n"
//						"		float bBD = By - kBD * Bx;\n"
//
//						"		float kEC = (Ey - Cy) * 1.0 / (Ex - Cx);\n"
//						"		float bEC = Ey - kEC * Ex;\n"
//
//						"		int flag = 0;\n"
//						"		float x = gl_TexCoord[0].x;\n"
//						"		float y = gl_TexCoord[0].y;\n"
//						"		if (y <= kAC * x + bAC && y <= kDA * x + bDA && (y >= kBD * x + bBD || y >= kEC * x + bEC))\n"
//						"		{\n"
//						"			flag = 1;\n"
//						"		}\n"
//						"		else if (y <= kEB * x + bEB)\n"
//						"		{\n"
//						"			if ((y >= kEC * x + bEC && y >= kDA * x + bDA) || (y >= kAC * x + bAC && y >= kBD * x + bBD))\n"
//						"			{\n"
//						"				flag = 1;\n"
//						"			}\n"
//						"			else\n"
//						"			{\n"
//						"				flag = 0;\n"
//						"			}\n"
//						"		}\n"
//						"		else\n"
//						"		{\n"
//						"			flag = 0;\n"
//						"		}\n"
//						"		if(flag<1)\n"
//						"			discard;\n"
//						"	}\n"
//						"	gl_FragColor = gl_Color;\n"
//						"}\n"
//						;
//					vertShader = 
//						"#version 330 compatibility\n"
//						"attribute vec4 a_color;\n"
//						"attribute float a_size;\n"
//						"attribute float a_rotation;\n"
//						"attribute float a_markerType;\n"
//
//						"uint oe_index_objectid; \n"                        // Vertex stage global containing the Object ID; set in ObjectIndex shader.
//						"flat out vec4 oe_pick_encoded_objectid; \n"        // output encoded oid to fragment shader
//						"flat out int  oe_pick_color_contains_objectid; \n" // whether color already contains oid (written by another RTT camera)
//
//						//"varying float v_size;\n"
//						//"varying float v_rotation;\n"
//						"void main(){\n"
//						//"	v_size = a_size;\n"
//						//"	v_rotation = a_rotation;\n"
//						"    oe_pick_color_contains_objectid = (oe_index_objectid == 1u) ? 1 : 0; \n"
//						"    if ( oe_pick_color_contains_objectid == 0 ) \n"
//						"    { \n"
//						"        float b0 = float((oe_index_objectid & 0xff000000u) >> 24u); \n"
//						"        float b1 = float((oe_index_objectid & 0x00ff0000u) >> 16u); \n"
//						"        float b2 = float((oe_index_objectid & 0x0000ff00u) >> 8u ); \n"
//						"        float b3 = float((oe_index_objectid & 0x000000ffu)       ); \n"
//						"        oe_pick_encoded_objectid = vec4(b0, b1, b2, b3) * 0.00392156862; \n" // i.e. 1/2558
//						"    } \n"
//						"	gl_TexCoord[0] = vec4(a_size,a_rotation,a_markerType, 1.0);\n"
//						"	gl_FrontColor = a_color;\n"
//						"	gl_Position = gl_Vertex;\n"
//						"}\n"
//						;
//					geoShader = 
//						"#version 330 compatibility\n"
//						"#extension GL_EXT_geometry_shader4:enable\n"
//						"layout(points) in;\n"
//						"layout(triangle_strip,max_vertices=1024) out;\n"
//
//						"uniform vec3 earthHeart;"
//						"void main(){\n"
//						//"	testVec = gl_ModelViewMatrix * gl_ProjectionMatrix * testVec;\n"
//						"	float r = sqrt(earthHeart.x*earthHeart.x+earthHeart.y*earthHeart.y+earthHeart.z*earthHeart.z);\n"
//						"	float a = gl_PositionIn[0].x-earthHeart.x;\n"
//						"	float b = gl_PositionIn[0].y-earthHeart.y;\n"
//						"	float c = gl_PositionIn[0].z-earthHeart.z;\n"
//						"	float w = gl_PositionIn[0].w;\n"
//
//
//						"	float s = gl_TexCoordIn[0][0].x/30000;\n"
//						"	float ra = radians(gl_TexCoordIn[0][0].y-45);\n"
//						"	float rb = radians(gl_TexCoordIn[0][0].y+45);\n"
//						"	float rc = radians(gl_TexCoordIn[0][0].y+135);\n"
//						"	float rd = radians(gl_TexCoordIn[0][0].y+225);\n"
//						"	float markerType = gl_TexCoordIn[0][0].z;\n"
//
//						"	if((r*r)/b-b>0)\n"
//						"	{\n"
//						"		gl_TexCoord[0] = vec4(1, 1, markerType, 1);\n "
//						"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rc)+r*c/b*sin(rc)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rc),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rc)-r*a/b*sin(rc)),w);\n"
//						"		gl_FrontColor = gl_FrontColorIn[0];\n"  
//						"		EmitVertex();\n" 
//
//						"		gl_TexCoord[0] = vec4(1, 0, markerType, 1);\n "
//						"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rb)+r*c/b*sin(rb)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rb),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rb)-r*a/b*sin(rb)),w);\n;"
//						"		gl_FrontColor = gl_FrontColorIn[0];\n"  
//						"		EmitVertex();\n"
//
//						"		gl_TexCoord[0] = vec4(0, 1, markerType, 1);\n "
//						"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rd)+r*c/b*sin(rd)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rd),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rd)-r*a/b*sin(rd)),w);\n"
//						"		gl_FrontColor = gl_FrontColorIn[0];\n"  
//						"		EmitVertex();\n"
//
//						"		gl_TexCoord[0] = vec4(0, 0, markerType, 1);\n "
//						"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(ra)+r*c/b*sin(ra)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(ra),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(ra)-r*a/b*sin(ra)),w);\n"
//						"		gl_FrontColor = gl_FrontColorIn[0];\n"  
//						"		EmitVertex();\n"
//						"	}\n"
//
//						"	else\n"
//						"	{\n"
//						"		gl_TexCoord[0] = vec4(1, 0, markerType, 1);\n "
//						"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rc)+r*c/b*sin(rc)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rc),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rc)-r*a/b*sin(rc)),w);\n"
//						"		gl_FrontColor = gl_FrontColorIn[0];\n"  
//						"		EmitVertex();\n"
//
//						"		gl_TexCoord[0] = vec4(1, 1, markerType, 1);\n "
//						"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rb)+r*c/b*sin(rb)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rb),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rb)-r*a/b*sin(rb)),w);\n;"
//						"		gl_FrontColor = gl_FrontColorIn[0];\n"  
//						"		EmitVertex();\n"
//
//						"		gl_TexCoord[0] = vec4(0, 0, markerType, 1);\n "
//						"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(rd)+r*c/b*sin(rd)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(rd),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(rd)-r*a/b*sin(rd)),w);\n"
//						"		gl_FrontColor = gl_FrontColorIn[0];\n"  
//						"		EmitVertex();\n"
//
//						"		gl_TexCoord[0] = vec4(0, 1, markerType, 1);\n "
//						"		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_PositionIn[0].x+s/(2*sqrt((r*r)/(b*b)-1))*(-a*cos(ra)+r*c/b*sin(ra)),gl_PositionIn[0].y+s/(2*sqrt((r*r)/(b*b)-1))*((r*r)/b-b)*cos(ra),gl_PositionIn[0].z+s/(2*sqrt((r*r)/(b*b)-1))*(-c*cos(ra)-r*a/b*sin(ra)),w);\n"
//						"		gl_FrontColor = gl_FrontColorIn[0];\n"  
//						"		EmitVertex();\n"
//						"	}\n"
//						"	EndPrimitive();\n"
//						"}\n"
//						;
//				program->addShader(new osg::Shader(osg::Shader::FRAGMENT,fragShader));
//				program->addShader(new osg::Shader(osg::Shader::VERTEX,vertShader));
//				program->addShader(new osg::Shader(osg::Shader::GEOMETRY,geoShader));
//				ss->setAttributeAndModes(program,osg::StateAttribute::ON);
//				osg::StateAttribute::GLModeValue disable = osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE | osg::StateAttribute::PROTECTED;
//
//				ss->setMode(GL_BLEND,     disable );    
//				ss->setMode(GL_LIGHTING,  disable );
//				ss->setMode(GL_CULL_FACE, disable );
//
//				// Disabling GL_BLEND is not enough, because osg::Text re-enables it
//				// without regard for the OVERRIDE.
//				ss->setAttributeAndModes(new osg::BlendFunc(GL_ONE, GL_ZERO), osg::StateAttribute::OVERRIDE);
//
//				// install the picking shaders:
//				/*VirtualProgram* vp = createRTTProgram();
//				ss->setAttribute( vp );*/
//
//				// designate this as a pick camera, overriding any defaults below
//				ss->addUniform( new osg::Uniform("oe_isPickCamera", true), osg::StateAttribute::OVERRIDE );
//
//				// default value for the objectid override uniform:
//				std::cout<<Registry::objectIndex()->getObjectIDUniformName()<<std::endl;
//				ss->addUniform( new osg::Uniform(Registry::objectIndex()->getObjectIDUniformName().c_str(), 0u) );
//
//				geode->addDrawable( osgGeom );
//				if ( context.featureIndex() )
//					context.featureIndex()->tagDrawable( osgGeom, input );
//
//				// install clamping attributes if necessary
//				if (_style.has<AltitudeSymbol>() &&
//					_style.get<AltitudeSymbol>()->technique() == AltitudeSymbol::TECHNIQUE_GPU)
//				{            
//					Clamping::applyDefaultClampingAttrs( osgGeom, input->getDouble("__oe_verticalOffset", 0.0) );
//				}
//			}
//		}
//	}
#pragma endregion   
    return geode;
}

#define CROP_POLYS_BEFORE_TESSELLATING 1

void
BuildGeometryFilter::tileAndBuildPolygon(Geometry*               ring,
                                         const SpatialReference* featureSRS,
                                         const SpatialReference* mapSRS,
                                         bool                    makeECEF,
                                         bool                    tessellate,
                                         osg::Geometry*          osgGeom,
                                         const osg::Matrixd      &world2local)
{
#ifdef CROP_POLYS_BEFORE_TESSELLATING

#define MAX_POINTS_PER_CROP_TILE 1024

    bool built = false;
    unsigned count = ring->getTotalPointCount();
    if ( count > MAX_POINTS_PER_CROP_TILE )
    {
        unsigned tiles = (count / MAX_POINTS_PER_CROP_TILE) + 1u;
        double tx = ceil(sqrt((double)tiles));
        double ty = tx;
        Bounds b = ring->getBounds();
        double tw = b.width() / tx;
        double th = b.height() / ty;

        OE_DEBUG << "Found " << count << " points; cropping to " << tx << " x " << ty << std::endl;

        // Get the average Z, since GEOS will set teh Z of new verts to that of the cropping polygon,
        // which is stupid but that's how it is.
        double z = 0.0;
        for(unsigned i=0; i<ring->size(); ++i)
            z += ring->at(i).z();
        z /= ring->size();

        osg::ref_ptr<Polygon> poly = new Polygon;
        poly->resize( 4 );

        built = true;
        for(int x=0; x<(int)tx; ++x)
        {
            for(int y=0; y<(int)ty; ++y)
            {
                (*poly)[0].set( b.xMin() + tw*(double)x,     b.yMin() + th*(double)y,     z );
                (*poly)[1].set( b.xMin() + tw*(double)(x+1), b.yMin() + th*(double)y,     z );
                (*poly)[2].set( b.xMin() + tw*(double)(x+1), b.yMin() + th*(double)(y+1), z );
                (*poly)[3].set( b.xMin() + tw*(double)x,     b.yMin() + th*(double)(y+1), z );
                
                osg::ref_ptr<Geometry> ringTile;
                if ( ring->crop(poly.get(), ringTile) )
                {
                    // Use an iterator since crop could return a multi-polygon
                    GeometryIterator gi( ringTile.get(), false );
                    while( gi.hasMore() )
                    {
                        Geometry* geom = gi.next();
                        buildPolygon(geom, featureSRS, mapSRS, makeECEF, tessellate, osgGeom, world2local);
                    }
                }
                else 
                {
                    // If crop resulted in empty geometry (valid case) ringTile will still be valid,
                    // otherwise we need to process the entire polygon without tiling.
                    if (!ringTile.valid())
                    {
                        //clean up geometry
                        osgGeom->setVertexArray(0L);
                        if (osgGeom->getNumPrimitiveSets())
                            osgGeom->removePrimitiveSet(0, osgGeom->getNumPrimitiveSets());

                        OE_NOTICE << LC << "GEOS crop failed, tessellating feature without tiling." << std::endl;

                        built = false;
                        break;
                    }
                }
            }

            // GEOS failed 
            if (!built)
                break;
        }
    }

    if ( !built )
    {
        buildPolygon(ring, featureSRS, mapSRS, makeECEF, tessellate, osgGeom, world2local);
    }
    

    if ( tessellate )
    {
        osgEarth::Tessellator oeTess;
        if (!oeTess.tessellateGeometry(*osgGeom))
        {
            //fallback to osg tessellator
            OE_DEBUG << LC << "Falling back on OSG tessellator (" << osgGeom->getName() << ")" << std::endl;

            osgUtil::Tessellator tess;
            tess.setTessellationType( osgUtil::Tessellator::TESS_TYPE_GEOMETRY );
            tess.setWindingType( osgUtil::Tessellator::TESS_WINDING_POSITIVE );
            tess.retessellatePolygons( *osgGeom );
        }
    }    
    
    osgUtil::SmoothingVisitor::smooth( *osgGeom );

#else

    // non-cropped way
    buildPolygon(ring, featureSRS, mapSRS, makeECEF, tessellate, osgGeom, world2local);
    if ( tessellate )
    {
        osgEarth::Tessellator oeTess;
        if (!oeTess.tessellateGeometry(*osgGeom))
        {
            //fallback to osg tessellator
            OE_INFO << LC << "OE Tessellation failed! Using OSG tessellator. (" << osgGeom->getName() << ")" << std::endl;

            osgUtil::Tessellator tess;
            tess.setTessellationType( osgUtil::Tessellator::TESS_TYPE_GEOMETRY );
            tess.setWindingType( osgUtil::Tessellator::TESS_WINDING_POSITIVE );
            tess.retessellatePolygons( *osgGeom );
        }
    }
#endif
}

// builds and tessellates a polygon (with or without holes)
void
BuildGeometryFilter::buildPolygon(Geometry*               ring,
                                  const SpatialReference* featureSRS,
                                  const SpatialReference* mapSRS,
                                  bool                    makeECEF,
                                  bool                    tessellate,
                                  osg::Geometry*          osgGeom,
                                  const osg::Matrixd      &world2local)
{
    if ( !ring->isValid() )
        return;

    ring->rewind(osgEarth::Symbology::Geometry::ORIENTATION_CCW);

    osg::ref_ptr<osg::Vec3Array> allPoints = new osg::Vec3Array();
    transformAndLocalize( ring->asVector(), featureSRS, allPoints.get(), mapSRS, world2local, makeECEF );

    Polygon* poly = dynamic_cast<Polygon*>(ring);
    if ( poly )
    {
        RingCollection ordered(poly->getHoles().begin(), poly->getHoles().end());
        std::sort(ordered.begin(), ordered.end(), holeCompare);

        for( RingCollection::const_iterator h = ordered.begin(); h != ordered.end(); ++h )
        {
            Geometry* hole = h->get();
            if ( hole->isValid() )
            {
                hole->rewind(osgEarth::Symbology::Geometry::ORIENTATION_CW);

                osg::ref_ptr<osg::Vec3Array> holePoints = new osg::Vec3Array();
                transformAndLocalize( hole->asVector(), featureSRS, holePoints.get(), mapSRS, world2local, makeECEF );

                // find the point with the highest x value
                unsigned int hCursor = 0;
                for (unsigned int i=1; i < holePoints->size(); i++)
                {
                    if ((*holePoints)[i].x() > (*holePoints)[hCursor].x())
                      hCursor = i;
                }

                double x1 = (*holePoints)[hCursor].x();
                double y1 = (*holePoints)[hCursor].y();
                double y2 = (*holePoints)[hCursor].y();

                unsigned int edgeCursor = UINT_MAX;
                double edgeDistance = DBL_MAX;
                unsigned int foundPointCursor = UINT_MAX;
                for (unsigned int i=0; i < allPoints->size(); i++)
                {
                    unsigned int next = i == allPoints->size() - 1 ? 0 : i + 1;
                    double xMax = osg::maximum((*allPoints)[i].x(), (*allPoints)[next].x());

                    if (xMax > (*holePoints)[hCursor].x())
                    {
                        double x2 = xMax + 1.0;
                        double x3 = (*allPoints)[i].x();
                        double y3 = (*allPoints)[i].y();
                        double x4 = (*allPoints)[next].x();
                        double y4 = (*allPoints)[next].y();

                        double xi=0.0, yi=0.0;
                        bool intersects = false;
                        unsigned int hitPointCursor = UINT_MAX;
                        if (y1 == y3 && x3 > x1)
                        {
                            xi = x3;
                            hitPointCursor = i;
                            intersects = true;
                        }
                        else if (y1 == y4 && x4 > x1)
                        {
                            xi = x4;
                            hitPointCursor = next;
                            intersects = true;
                        }
                        else if (segmentsIntersect(x1, y1, x2, y2, x3, y3, x4, y4, xi, yi))
                        {
                            intersects = true;
                        }

                        double dist = (osg::Vec2d(xi, yi) - osg::Vec2d(x1, y1)).length();
                        if (intersects && dist < edgeDistance)
                        {
                            foundPointCursor = hitPointCursor;
                            edgeCursor = hitPointCursor != UINT_MAX ? hitPointCursor : (x3 >= x4 ? i : next);
                            edgeDistance = dist;
                        }
                    }
                }

                if (foundPointCursor == UINT_MAX && edgeCursor != UINT_MAX)
                {
                    // test for intersecting edges between x1 and x2
                    // (skipping the two segments for which edgeCursor is a vert)

                    double x2 = (*allPoints)[edgeCursor].x();
                    y2 = (*allPoints)[edgeCursor].y();

                    bool foundIntersection = false;
                    for (unsigned int i=0; i < allPoints->size(); i++)
                    {
                        unsigned int next = i == allPoints->size() - 1 ? 0 : i + 1;

                        if (i == edgeCursor || next == edgeCursor)
                          continue;

                        double x3 = (*allPoints)[i].x();
                        double y3 = (*allPoints)[i].y();
                        double x4 = (*allPoints)[next].x();
                        double y4 = (*allPoints)[next].y();

                        foundIntersection = foundIntersection || segmentsIntersect(x1, y1, x2, y2, x3, y3, x4, y4);

                        if (foundIntersection)
                        {
                            unsigned int prev = i == 0 ? allPoints->size() - 1 : i - 1;

                            if (!isCCW((*allPoints)[prev].x(), (*allPoints)[prev].y(), x3, y3, x4, y4))
                            {
                                edgeCursor = i;
                                x2 = (*allPoints)[edgeCursor].x();
                                y2 = (*allPoints)[edgeCursor].y();
                                foundIntersection = false;
                            }
                        }

                    }
                }

                if (edgeCursor != UINT_MAX)
                {
                    // build array of correctly ordered new points to add to the outer loop
                    osg::ref_ptr<osg::Vec3Array> insertPoints = new osg::Vec3Array();
                    insertPoints->reserve(holePoints->size() + 2);

                    unsigned int p = hCursor;
                    do
                    {
                        insertPoints->push_back((*holePoints)[p]);
                        p = p == holePoints->size() - 1 ? 0 : p + 1;
                    } while(p != hCursor);

                    insertPoints->push_back((*holePoints)[hCursor]);
                    insertPoints->push_back((*allPoints)[edgeCursor]);
                    
                    // insert new points into outer loop
                    osg::Vec3Array::iterator it = edgeCursor == allPoints->size() - 1 ? allPoints->end() : allPoints->begin() + (edgeCursor + 1);
                    allPoints->insert(it, insertPoints->begin(), insertPoints->end());
                }
            }
        }
    }
    
    GLenum mode = GL_LINE_LOOP;
    if ( osgGeom->getVertexArray() == 0L )
    {
        osgGeom->addPrimitiveSet( new osg::DrawArrays( mode, 0, allPoints->size() ) );
        osgGeom->setVertexArray( allPoints.get() );
    }
    else
    {
        osg::Vec3Array* v = static_cast<osg::Vec3Array*>(osgGeom->getVertexArray());
        osgGeom->addPrimitiveSet( new osg::DrawArrays( mode, v->size(), allPoints->size() ) );
        //v->reserve(v->size() + allPoints->size());
        std::copy(allPoints->begin(), allPoints->end(), std::back_inserter(*v));
    }

    //// Normal computation.
    //// Not completely correct, but better than no normals at all. TODO: update this
    //// to generate a proper normal vector in ECEF mode.
    ////
    //// We cannot accurately rely on triangles from the tessellation, since we could have
    //// very "degraded" triangles (close to a line), and the normal computation would be bad.
    //// In this case, we would have to average the normal vector over each triangle of the polygon.
    //// The Newell's formula is simpler and more direct here.
    //osg::Vec3 normal( 0.0, 0.0, 0.0 );
    //for ( size_t i = 0; i < poly->size(); ++i )
    //{
    //    osg::Vec3 pi = (*poly)[i];
    //    osg::Vec3 pj = (*poly)[ (i+1) % poly->size() ];
    //    normal[0] += ( pi[1] - pj[1] ) * ( pi[2] + pj[2] );
    //    normal[1] += ( pi[2] - pj[2] ) * ( pi[0] + pj[0] );
    //    normal[2] += ( pi[0] - pj[0] ) * ( pi[1] + pj[1] );
    //}
    //normal.normalize();

    //osg::Vec3Array* normals = new osg::Vec3Array();
    //normals->push_back( normal );
    //osgGeom->setNormalArray( normals );
    //osgGeom->setNormalBinding( osg::Geometry::BIND_OVERALL );
}


osg::Node*
BuildGeometryFilter::push( FeatureList& input, FilterContext& context )
{
    osg::ref_ptr<osg::Group> result = new osg::Group();

    computeLocalizers( context );

    const LineSymbol*    line  = _style.get<LineSymbol>();
    const PolygonSymbol* poly  = _style.get<PolygonSymbol>();
    const PointSymbol*   point = _style.get<PointSymbol>();

    // bin the feautres into polygons, lines, polygonized lines, and points.
    FeatureList polygons;
    FeatureList lines;
    FeatureList polygonizedLines;
    FeatureList points;

    for(FeatureList::iterator i = input.begin(); i != input.end(); ++i)
    {
        Feature* f = i->get();

        // first consider the overall style:
        bool has_polysymbol     = poly != 0L;
        bool has_linesymbol     = line != 0L && line->stroke()->widthUnits() == Units::PIXELS;
        bool has_polylinesymbol = line != 0L && line->stroke()->widthUnits() != Units::PIXELS;
        bool has_pointsymbol    = point != 0L;

        // if the featue has a style set, that overrides:
        if ( f->style().isSet() )
        {
            has_polysymbol     = has_polysymbol     || (f->style()->has<PolygonSymbol>());
            has_linesymbol     = has_linesymbol     || (f->style()->has<LineSymbol>() && f->style()->get<LineSymbol>()->stroke()->widthUnits() == Units::PIXELS);
            has_polylinesymbol = has_polylinesymbol || (f->style()->has<LineSymbol>() && f->style()->get<LineSymbol>()->stroke()->widthUnits() != Units::PIXELS);
            has_pointsymbol    = has_pointsymbol    || (f->style()->has<PointSymbol>());
        }

        // if no style is set, use the geometry type:
        if ( !has_polysymbol && !has_linesymbol && !has_polylinesymbol && !has_pointsymbol && f->getGeometry() )
        {
            switch( f->getGeometry()->getComponentType() )
            {
            case Geometry::TYPE_LINESTRING:
            case Geometry::TYPE_RING:
                f->style()->add( new LineSymbol() );
                has_linesymbol = true;
                break;

            case Geometry::TYPE_POINTSET:
                f->style()->add( new PointSymbol() );
                has_pointsymbol = true;
                break;

            case Geometry::TYPE_POLYGON:
                f->style()->add( new PolygonSymbol() );
                has_polysymbol = true;
                break;
            }
        }

        if ( has_polysymbol )
            polygons.push_back( f );

        if ( has_linesymbol )
            lines.push_back( f );

        if ( has_polylinesymbol )
            polygonizedLines.push_back( f );

        if ( has_pointsymbol )
            points.push_back( f );
    }

    // process them separately.

    if ( polygons.size() > 0 )
    {
        OE_TEST << LC << "Building " << polygons.size() << " polygons." << std::endl;
        osg::ref_ptr<osg::Geode> geode = processPolygons(polygons, context);
        if ( geode->getNumDrawables() > 0 )
        {
            osgUtil::Optimizer o;
            o.optimize( geode.get(), 
                osgUtil::Optimizer::MERGE_GEOMETRY |
                osgUtil::Optimizer::INDEX_MESH |
                osgUtil::Optimizer::VERTEX_PRETRANSFORM |
                osgUtil::Optimizer::VERTEX_POSTTRANSFORM );

            result->addChild( geode.get() );
        }
    }

    if ( polygonizedLines.size() > 0 )
    {
        OE_TEST << LC << "Building " << polygonizedLines.size() << " polygonized lines." << std::endl;
        bool twosided = polygons.size() > 0 ? false : true;
        osg::ref_ptr<osg::Geode> geode = processPolygonizedLines(polygonizedLines, twosided, context);
        if ( geode->getNumDrawables() > 0 )
        {
            osgUtil::Optimizer o;
            o.optimize( geode.get(), 
                osgUtil::Optimizer::MERGE_GEOMETRY |
                osgUtil::Optimizer::INDEX_MESH |
                osgUtil::Optimizer::VERTEX_PRETRANSFORM |
                osgUtil::Optimizer::VERTEX_POSTTRANSFORM );

            result->addChild( geode.get() );
        }
    }

    if ( lines.size() > 0 )
    {
        OE_TEST << LC << "Building " << lines.size() << " lines." << std::endl;
        osg::ref_ptr<osg::Geode> geode = processLines(lines, context);
        if ( geode->getNumDrawables() > 0 )
        {
            osgUtil::Optimizer o;
            o.optimize( geode.get(), 
                osgUtil::Optimizer::MERGE_GEOMETRY );

            applyLineSymbology( geode->getOrCreateStateSet(), line );
            result->addChild( geode.get() );
        }
    }

    if ( points.size() > 0 )
    {
        OE_TEST << LC << "Building " << points.size() << " points." << std::endl;
        osg::ref_ptr<osg::Geode> geode = processPoints(points, context);
        if ( geode->getNumDrawables() > 0 )
        {
            osgUtil::Optimizer o;
            o.optimize( geode.get(), 
                osgUtil::Optimizer::MERGE_GEOMETRY );

            applyPointSymbology( geode->getOrCreateStateSet(), point );
            result->addChild( geode.get() );
        }
    }

    // indicate that geometry contains clamping attributes
    if (_style.has<AltitudeSymbol>() &&
        _style.get<AltitudeSymbol>()->technique() == AltitudeSymbol::TECHNIQUE_GPU)
    {
        Clamping::installHasAttrsUniform( result->getOrCreateStateSet() );
    }    

    // Prepare buffer objects.
    AllocateAndMergeBufferObjectsVisitor allocAndMerge;
    result->accept( allocAndMerge );


    if ( result->getNumChildren() > 0 )
    {
        // apply the delocalization matrix for no-jitter
        return delocalize( result.release() );
    }
    else
    {
        return 0L;
    }
}

osg::Image* BuildGeometryFilter::createImage(osg::Vec4Array* pointVec)
{
	osg::Image* image = new osg::Image;
	int numpixels = pointVec->size();
	image->allocateImage(numpixels,1,1,GL_RGBA,GL_FLOAT);
	image->setInternalTextureFormat(GL_RGBA);
	osg::Vec4f* datePtr = (osg::Vec4f*)image->data();
	for (int i=0;i<numpixels;i++)
	{
		*datePtr++ = pointVec->at(i);
	}
	return image;
}

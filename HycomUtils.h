/**********************************************************/
//name          : HycomUtils.h
//function      : HycomÁ÷³¡
//copyright     : 
//author        : JXF
//date          : 2017-6-1
/**********************************************************/
#ifndef HYCOMUTILS_H
#define HYCOMUTILS_H

#include <string>
#include <QDate>
#include <osgEarth/GeoData>
#include <osgEarthAnnotation/PlaceNode>
#include <osgEarthAnnotation/FeatureNode>
#include <osgEarthUtil/AnnotationEvents>
#include <osgEarthAnnotation/HighlightDecoration>
#include <osgEarthAnnotation/FeatureNode>
#include <QSqlDatabase>
#include <QWidget>
#include <osgEarthUtil/ObjectLocator>
#include <osg/MatrixTransform>

using namespace osgEarth;
using namespace osgEarth::Annotation;
using namespace osgEarth::Util;
using namespace std;

class HycompathInfoGroup : public osg::Group
{
private:
	//static HycompathInfoGroup* hycompathInfoGroup;
	string _name;
	string _platformId;

	struct PathInfo
	{
		double poingLength;
		double sumLength;
		string pointID;
		double lat;
		double lon;
		double speed;
	};

	struct LineInfo
	{
		osg::ref_ptr<FeatureNode> lineNode;
		vector<PathInfo> pathInfoVec;
		double allLength;
		double nowDistance;
		int nowIndex;
		bool done;
	};

	vector<LineInfo> _lineInfoVec;
	
	Style _argoPathLineStyle;
	void _initialize();

	osg::Group* _argoCycGroup;
public:
	static HycompathInfoGroup* getInstance();
	/*void replay();*/
	HycompathInfoGroup();
	~HycompathInfoGroup();
	/*void removeHycompathInfoGroup();*/
	/*virtual void traverse(osg::NodeVisitor& nv);*/
};



#endif
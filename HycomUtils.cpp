/**********************************************************/
//name          : HycomUtils.cpp
//function      : Hycom流场
//copyright     : 
//author        : JXF
//date          : 2017-6-1
/**********************************************************/

#include "HycomUtils.h"
#include "ReferenceTheInstances.h"
#include <osgUtil/Optimizer>
#include <osgEarth/GeoMath>
#include <QMessageBox>
#include "SystemSettingUtils.h"
#include "StringUtils.h"
#define D2R (osg::PI/180.0)
#define R2D (180.0/osg::PI)

//HycompathInfoGroup* HycompathInfoGroup::hycompathInfoGroup = NULL;

HycompathInfoGroup* HycompathInfoGroup::getInstance()
{
	//if (hycompathInfoGroup != NULL)
	//{
	//	//hycompathInfoGroup = new HycompathInfoGroup();
	//	//delete hycompathInfoGroup;
	//	hycompathInfoGroup=NULL;
	//}
	//hycompathInfoGroup = new HycompathInfoGroup();
	return new HycompathInfoGroup();
}

HycompathInfoGroup::HycompathInfoGroup()
{
	//ReferenceTheInstances::MapNodePtr->addChild(this);

	//_argoCycGroup = new osg::Group();
	//ReferenceTheInstances::MapNodePtr->addChild(_argoCycGroup);
	/*_argoPathLineStyle.getOrCreate<LineSymbol>()->stroke()->color() = Color::White;
	 _argoPathLineStyle.getOrCreate<LineSymbol>()->trackLineString() = "*";
	_argoPathLineStyle.getOrCreate<LineSymbol>()->stroke()->width() = 2.0f;
	_argoPathLineStyle.getOrCreate<AltitudeSymbol>()->clamping() = AltitudeSymbol::CLAMP_TO_TERRAIN;
	_argoPathLineStyle.getOrCreate<AltitudeSymbol>()->technique() = AltitudeSymbol::TECHNIQUE_DRAPE;*/
	_argoPathLineStyle.getOrCreate<PointSymbol>()->isMutliple() = true;
	_argoPathLineStyle.getOrCreate<PointSymbol>()->isFlow() = true;
	//_argoPathLineStyle.getOrCreate<PointSymbol>()->isProfile() = true;
	_argoPathLineStyle.getOrCreate<AltitudeSymbol>()->clamping() = AltitudeSymbol::CLAMP_TO_TERRAIN;
	_argoPathLineStyle.getOrCreate<AltitudeSymbol>()->technique() = AltitudeSymbol::TECHNIQUE_GPU;

	_initialize();
}

HycompathInfoGroup::~HycompathInfoGroup()
{
	
}


//void HycompathInfoGroup::removeHycompathInfoGroup()
//{
//	_done = true;
//	this->removeChildren(0,this->getNumChildren());
//	_argoCycGroup->removeChildren(0,_argoCycGroup->getNumChildren());
//}

//void HycompathInfoGroup::replay()
//{
//	if (_pathInfoVec.size() == 0)
//	{
//		int ret = QMessageBox::warning(0, QString::fromLocal8Bit("警告"),QString::fromLocal8Bit("轨迹数据缺失"));
//		this->removeHycompathInfoGroup();
//		return;
//	}
//	double originLon = _pathInfoVec[0].pointNode->getFeature()->getGeometry()->at(0).x();
//	double originLat = _pathInfoVec[0].pointNode->getFeature()->getGeometry()->at(0).y();
//	_originNode->setPosition(GeoPoint(ReferenceTheInstances::MapNodePtr->getMapSRS()->getGeographicSRS(),originLon,originLat));
//	_initialize();
//}

//void HycompathInfoGroup::_initialize()
//{
//	_lineInfoVec.clear();
//	this->removeChildren(0,this->getNumChildren());
//	_argoCycGroup->removeChildren(0,_argoCycGroup->getNumChildren());
//	/*string sql = "select distinct lineid,speed from hycom.hycom_glb_912_2016050100_t000_uv3z_0_pppp limit 150";
//	QSqlQuery query(QString::fromStdString(sql),ReferenceTheInstances::db);*/
//	QSqlQuery query1(ReferenceTheInstances::db);
//	
//	string lastLineId = "";
//	string sql = "select lineid,pointid,lat,lon,speed from hycom.hycom_glb_912_2016050100_t000_uv3z_0_pppp where lineid<25 order by lineid,pointid";
//	query1.exec(QString::fromStdString(sql));
//	int count = 0;
//	LineInfo* currentLineInfo;
//	LineInfo lineInfo;
//	while(query1.next())
//	{
//		if (lastLineId!=query1.value(0).toString().toStdString())
//		{		
//			if (lastLineId!="")
//			{
//				_lineInfoVec.push_back(*currentLineInfo);
//			}
//			lineInfo.nowDistance = 0.0;
//			lineInfo.nowIndex = 0;
//			lineInfo.done = false;
//			lineInfo.allLength = 0;
//			lineInfo.pathInfoVec.clear();
//			/*osgEarth::Symbology::LineString* pathLine = new osgEarth::Symbology::LineString();
//			pathLine->push_back(query1.value(3).toDouble(),query1.value(2).toDouble());
//			pathLine->push_back(120,30);*/
//
//			osg::ref_ptr<osgEarth::Symbology::PointSet> pointSet = new osgEarth::Symbology::PointSet();
//			pointSet->push_back(query1.value(3).toDouble(),query1.value(2).toDouble());
//			lineInfo.lineNode = new FeatureNode(ReferenceTheInstances::MapNodePtr, new Feature(pointSet, ReferenceTheInstances::MapNodePtr->getMapSRS()->getGeographicSRS(), _argoPathLineStyle));
//			this->addChild( lineInfo.lineNode );
//			count = 0;
//			lastLineId = query1.value(0).toString().toStdString();
//			currentLineInfo = &lineInfo;
//		}
//
//		PathInfo info;
//		info.pointID = query1.value(1).toDouble();
//		info.lat = query1.value(2).toDouble();
//		info.lon = query1.value(3).toDouble();
//		info.speed = query1.value(4).toDouble();
//
//		if(count==0)
//		{
//			info.poingLength = 0;
//			info.sumLength = 0;
//		}
//		else
//		{
//			double length = GeoMath::distance(info.lat,info.lon,currentLineInfo->pathInfoVec[count-1].lat,currentLineInfo->pathInfoVec[count-1].lon);
//			currentLineInfo->allLength += length;
//			info.poingLength = length;
//			info.sumLength = length + currentLineInfo->pathInfoVec[count-1].sumLength;
//		}
//		currentLineInfo->pathInfoVec.push_back(info);
//		count++;
//	}
//	if (count!=0)
//	{
//		_lineInfoVec.push_back(*currentLineInfo);
//	}
//	
//	if (_lineInfoVec.size() == 0)
//	{
//		int ret = QMessageBox::warning(0, QString::fromLocal8Bit("警告"),QString::fromLocal8Bit("轨迹数据缺失"));
//		return;
//
//	}
//}

void HycompathInfoGroup::_initialize()
{
	_lineInfoVec.clear();
	this->removeChildren(0,this->getNumChildren());
	
	string lastLineId = "";
	string sizeString = "*";
	string rotationString = "*";
	double count = 0;
	int totalCount = 0;
	double allcount = 0;
	//LineInfo* currentLineInfo;
	LineInfo lineInfo;
	lineInfo.lineNode=NULL;
	osg::ref_ptr<osgEarth::Symbology::PointSet> pointSet = new osgEarth::Symbology::PointSet();

	string hycompath = SystemSettingUtils::getMapSettingByNodePath("/mapSetting/hycomcsvpath");
	QFile file(QString::fromLocal8Bit(hycompath.c_str()));
	if(!file.open(QIODevice::ReadOnly))
		qDebug()<<"OPEN FILE FAILED";
	QTextStream* out = new QTextStream(&file);
	QString lineStr = out->readLine();
	//QStringList tempOption = out->readAll().split("\n");
	QStringList everyLine;
	
	while(lineStr!="")
	{
		everyLine = QString::QString(lineStr).remove("\r").split(",");
		if (lastLineId!=everyLine.at(0).toStdString())
		{	
			if (allcount>2361000)
			{
				_argoPathLineStyle.getOrCreate<PointSymbol>()->isFeatureNode() = true;
				_argoPathLineStyle.getOrCreate<PointSymbol>()->featureCount() = totalCount;
				_argoPathLineStyle.getOrCreate<PointSymbol>()->sizeString() = sizeString;
				_argoPathLineStyle.getOrCreate<PointSymbol>()->rotationString() = rotationString;
				lineInfo.lineNode = new FeatureNode(ReferenceTheInstances::MapNodePtr, new Feature(pointSet, ReferenceTheInstances::MapNodePtr->getMapSRS()->getGeographicSRS(), _argoPathLineStyle));
				this->addChild( lineInfo.lineNode );

				totalCount = 0;
				sizeString = "*";
				rotationString = "*";
				allcount = 0;
				//pointSet = new osgEarth::Symbology::PointSet();
				pointSet->clear();
			}

			if (allcount!=0)
			{
				totalCount++;
				sizeString += StringUtils::getStringFromDouble(count,0);sizeString += "*";
			}
			lastLineId = everyLine.at(0).toStdString();
			count = 0;
		}
		allcount++;
		pointSet->push_back(everyLine.at(2).toDouble(),everyLine.at(1).toDouble());
		rotationString += StringUtils::getStringFromDouble(everyLine.at(3).toDouble(),1);rotationString += "*";
		count++;
		lineStr = out->readLine();
	}
	_argoPathLineStyle.getOrCreate<PointSymbol>()->isFeatureNode() = true;
	_argoPathLineStyle.getOrCreate<PointSymbol>()->featureCount() = totalCount;
	_argoPathLineStyle.getOrCreate<PointSymbol>()->sizeString() = sizeString;
	_argoPathLineStyle.getOrCreate<PointSymbol>()->rotationString() = rotationString;
	lineInfo.lineNode = new FeatureNode(ReferenceTheInstances::MapNodePtr, new Feature(pointSet, ReferenceTheInstances::MapNodePtr->getMapSRS()->getGeographicSRS(), _argoPathLineStyle));
	this->addChild( lineInfo.lineNode );
	file.close();

	//string sql = "select lineid,lat,lon,speed from hycom.hycom_glb_912_2016050100_t000_uv3z_0_p11 order by lineid,pointid";
	//query1.exec(QString::fromStdString(sql));
	//while(query1.next())
	//{
	//	
	//	if (lastLineId!=query1.value(0).toString().toStdString())
	//	{	
	//		if (allcount>2361000)
	//		{
	//			_argoPathLineStyle.getOrCreate<PointSymbol>()->isFeatureNode() = true;
	//			_argoPathLineStyle.getOrCreate<PointSymbol>()->featureCount() = totalCount;
	//			_argoPathLineStyle.getOrCreate<PointSymbol>()->sizeString() = sizeString;
	//			_argoPathLineStyle.getOrCreate<PointSymbol>()->rotationString() = rotationString;
	//			lineInfo.lineNode = new FeatureNode(ReferenceTheInstances::MapNodePtr, new Feature(pointSet, ReferenceTheInstances::MapNodePtr->getMapSRS()->getGeographicSRS(), _argoPathLineStyle));
	//			this->addChild( lineInfo.lineNode );

	//			totalCount = 0;
	//			sizeString = "*";
	//			rotationString = "*";
	//			allcount = 0;
	//			//pointSet = new osgEarth::Symbology::PointSet();
	//			pointSet->clear();
	//		}
	//		
	//		if (allcount!=0)
	//		{
	//			totalCount++;
	//			sizeString += StringUtils::getStringFromDouble(count,0);sizeString += "*";
	//		}
	//		lastLineId = query1.value(0).toString().toStdString();
	//		count = 0;
	//	}
	//	allcount++;
	//	pointSet->push_back(query1.value(2).toDouble(),query1.value(1).toDouble());
	//	rotationString += StringUtils::getStringFromDouble(query1.value(3).toDouble(),1);rotationString += "*";
	//	count++;
	//}
	//
	//_argoPathLineStyle.getOrCreate<PointSymbol>()->isFeatureNode() = true;
	//_argoPathLineStyle.getOrCreate<PointSymbol>()->featureCount() = totalCount;
	//_argoPathLineStyle.getOrCreate<PointSymbol>()->sizeString() = sizeString;
	//_argoPathLineStyle.getOrCreate<PointSymbol>()->rotationString() = rotationString;
	//lineInfo.lineNode = new FeatureNode(ReferenceTheInstances::MapNodePtr, new Feature(pointSet, ReferenceTheInstances::MapNodePtr->getMapSRS()->getGeographicSRS(), _argoPathLineStyle));
	//this->addChild( lineInfo.lineNode );
}

//void HycompathInfoGroup::traverse(osg::NodeVisitor& nv)
//{
//	for (int i=0;i<_lineInfoVec.size();i++)
//	{
//		LineInfo* lineInfo = &_lineInfoVec[i];
//		if(!lineInfo->done)
//		{
//			double increaseRate = abs(lineInfo->pathInfoVec[lineInfo->nowIndex].speed*5000000);
//			lineInfo->nowDistance += increaseRate;
//
//			double lengthRate;
//			LineString* lineString = (LineString*)lineInfo->lineNode->getFeature()->getGeometry();
//			if(lineString->size()>0)
//			{
//				lineString->pop_back();
//			}
//
//			if(lineInfo->nowDistance>lineInfo->pathInfoVec[lineInfo->nowIndex].sumLength)
//			{
//				lineString->push_back(lineInfo->pathInfoVec[lineInfo->nowIndex].lon,
//					lineInfo->pathInfoVec[lineInfo->nowIndex].lat);
//
//				if(lineInfo->nowIndex>= lineInfo->pathInfoVec.size()-1)
//					lineInfo->done = true;
//				lineInfo->nowIndex++;
//			}
//			if(!lineInfo->done)
//			{
//				double outlon,outlat;
//				double nowLon = lineInfo->pathInfoVec[lineInfo->nowIndex].lon;
//				double nowLat = lineInfo->pathInfoVec[lineInfo->nowIndex].lat;
//				double oldLon = lineInfo->pathInfoVec[lineInfo->nowIndex-1].lon;
//				double oldLat = lineInfo->pathInfoVec[lineInfo->nowIndex-1].lat;
//				lengthRate = (lineInfo->nowDistance-lineInfo->pathInfoVec[lineInfo->nowIndex-1].sumLength)/lineInfo->pathInfoVec[lineInfo->nowIndex].poingLength;
//				GeoMath::interpolate(D2R*oldLat,D2R*oldLon,
//					D2R*nowLat,D2R*nowLon,
//					lengthRate,outlat,outlon);
//				if (lengthRate>1)
//				{
//					lineInfo->nowDistance = lineInfo->pathInfoVec[lineInfo->nowIndex].sumLength;
//					lineString->push_back(nowLon,nowLat);
//					lineInfo->lineNode->getFeature()->setGeometry(lineString);
//					lineInfo->lineNode->dirty();
//				}
//				else
//				{
//					lineString->push_back(R2D*outlon,R2D*outlat);
//					lineInfo->lineNode->getFeature()->setGeometry(lineString);
//					lineInfo->lineNode->dirty();
//				}
//			}
//			else
//			{
//				LineString* lineString1 = new LineString();
//				lineInfo->lineNode->getFeature()->setGeometry(lineString1);
//				lineInfo->lineNode->dirty();
//				lineInfo->nowDistance = 0.0;
//				lineInfo->nowIndex = 0;
//				lineInfo->done = false;
//				lineInfo->allLength = 0;
//			}
//		}
//	}
//	Group::traverse(nv);
//}
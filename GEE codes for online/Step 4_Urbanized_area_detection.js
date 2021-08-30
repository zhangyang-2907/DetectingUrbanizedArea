//----------------------------------------------------------------------------------------
//                      Step 4: Detecting urbanized areas
//----------------------------------------------------------------------------------------
//Note: there are some places that need to set manually, indicated by (******need to set manually******) 
//in Lines 8-9, 12-13, 274-276, 362, 427, 463, 491, and 505. 

//duration : 2000-2019 (******need to set manually******)
var startYear = 2000;
var endYear = 2019;

//Temporal period of NDVI : May to October (******need to set manually******) 
var startMonth_NDVI = 5;
var durationMonth_NDVI = 6;

//________________________________________________________________________________________
//                   study area (Taken Fuzhou city as an example)
//________________________________________________________________________________________
var roi = /* color: #fcfafa */ee.Geometry.Polygon(
        [[[119.09, 26.18],
          [119.53, 26.18],
          [119.53, 25.845],
          [119.09, 25.845]]]);   
          
Map.centerObject(roi, 10);

//________________________________________________________________________________________
//                                      functions
//________________________________________________________________________________________
//(1) Masking cloud pixels in Landsat SR data
var maskClouds = function(image) {
                var qa = image.select('pixel_qa');   
                var mask = qa.bitwiseAnd(2).eq(2).or(qa.bitwiseAnd(4).eq(4));
                var clean = image.updateMask(mask);
                return clean;
                };

var dataPreprocessing = function (dataName,startDate,endDate){
    var landsat_collection = ee.ImageCollection(dataName);
    var landsat_studyArea = landsat_collection.filterBounds(roi);
    var landsat = landsat_studyArea.filterDate(startDate,endDate);
    var landsat_imageNoClouds = landsat.map(maskClouds); 
    return landsat_imageNoClouds.map(function(img){return img.clip(roi);});
};

var indicators = function (dataSet,BLUE,GREEN,RED,NIR,SWIR1,SWIR2){
   dataSet = dataSet.map(function(image) {
             var ndvi = image.normalizedDifference([NIR, RED]);
             return  image.addBands(ndvi.rename("NDVI")).toFloat();
          });
        
   dataSet = dataSet.map(function(image) {
             var awei_sh = image.expression('((BLUE+2.5*GREEN-1.5*(NIR+SWIR1)-0.25*SWIR2)/10000)',
             {
               'NIR': image.select(NIR),
               'BLUE': image.select(BLUE),
               'SWIR1': image.select(SWIR1),
               'SWIR2': image.select(SWIR2),
               'GREEN': image.select(GREEN)
             });                    
             return  image.addBands(awei_sh.rename("AWEI_sh")).toFloat();  
          });  
     
   dataSet = dataSet.map(function(image) {
             var swir = image.expression('(SWIR2/10000)',
             {
                'SWIR2': image.select(SWIR2)
              });
             return  image.addBands(swir.rename("SWIR")).toFloat();  
          });
          
    return dataSet.select('NDVI','AWEI_sh','SWIR');
};

var annualMean = function (dataSet,indicator,startYear,endYear,startMonth,durationMonth){
    var years = ee.List.sequence(startYear,endYear);  
    var meanAnnual_indi = years.map(function(y) {
        var start = ee.Date.fromYMD(y, startMonth, 1);  
        var stop = start.advance(durationMonth, 'month');  
        var collection = dataSet.select([indicator]).filterDate(start, stop); 
        var mean = collection.reduce(ee.Reducer.mean()).set('year', y);
        return mean.addBands(mean.metadata('year'));
    });
    return ee.ImageCollection(meanAnnual_indi);
    
};

var Trend = function(dataSet,indicator){
  //1. senslope calculation
  var senslope = dataSet.reduce(ee.Reducer.sensSlope());
  
  //2. Mann-Kendall trend test
  
  //the valid pixel numbers for each pixel(excluding the masked pixels)
  var N_valid = dataSet.select(indicator).count();
  
  var ones = ee.Image.constant(1).clip(roi);//  create ones image   
  var zeros = ee.Image.constant(0).clip(roi); //  create zeros image
  
  dataSet = dataSet.toList(dataSet.size());//converting collection to list
  var N = dataSet.length();//image number of dataset
  
  var countArray = ee.List.sequence(0,N.subtract(1))
    .map(function (xi) {return ee.Image.constant(-1).clip(roi).set('countArray_xi',xi);});
  
  countArray = countArray.add(zeros);
  countArray = countArray.add(zeros);
  countArray = countArray.add(zeros);
  
 
  var first = countArray.add(ones);
  first = first.add(zeros);
 
  var dataSet_xi = dataSet.slice(0,N.subtract(1));
 
  var accumulateTies = function(image,list){
          
          var xi = ee.List(dataSet).indexOf(image);
          var inputArray_xi = image.select(indicator);
          var first_xj = ee.List(list).set(-1, inputArray_xi);
          
          var dataSet_xj = dataSet.slice(ee.Number(xi).add(1),dataSet.length());
          
          var accumulateTies_xj = function(image_xj,list_xj){
             var xj = ee.List(dataSet).indexOf(image_xj);
             var currentImage = image_xj.select(indicator);
             var previousImage = ee.Image(ee.List(list_xj).get(-1));
             var diff = (currentImage.subtract(previousImage));
            
             var countArray_xi = ee.Image(ee.List(list_xj).get(xi));
             var countArray_xj = ee.Image(ee.List(list_xj).get(xj));
             var count = ee.Image(ee.List(list_xj).get(-2));
             
             var ties = diff.eq(0).and(countArray_xi.neq(0));
             countArray_xj = countArray_xj.where(ties, zeros);
             count = count.where(ties,count.add(1));
             list_xj = ee.List(list_xj).set(xj,countArray_xj);
             list_xj = ee.List(list_xj).set(-2,count);
            
             var c_diff = diff.where(diff.gt(0),1).where(diff.neq(1),0);
             var concordant = ee.Image(ee.List(list_xj).get(-5));
             concordant = concordant.add(c_diff.unmask());
             list_xj = ee.List(list_xj).set(-5,concordant);
             
             var d_diff = diff.where(diff.eq(0),0).where(diff.gt(0),0).where(diff.lt(0),1);
             var discordant = ee.Image(ee.List(list_xj).get(-4));
             discordant = discordant.add(d_diff.unmask());
             list_xj = ee.List(list_xj).set(-4,discordant);
             
             var cd_diff = diff.where(diff.eq(0),0).where(diff.gt(0),1).where(diff.lt(0),-1);
             var cd = ee.Image(ee.List(list_xj).get(-3));
             cd = cd.add(cd_diff.unmask());
             list_xj = ee.List(list_xj).set(-3,cd);
    
             return ee.List(list_xj);
             
          };
          list = ee.ImageCollection(dataSet_xj).iterate(accumulateTies_xj,first_xj);
          
          list = ee.List(list).set(xi,ee.List(list).get(-2));
          list = ee.List(list).set(-2,ee.Image.constant(1).clip(roi));
          return ee.List(list);
  };
  var tieCount = ee.ImageCollection(dataSet_xi).iterate(accumulateTies,first);
  
  var tie_xi = ee.ImageCollection(ee.List(tieCount).slice(0,-6));
 
  var m0 = ee.Image(N_valid).multiply(ee.Image(N_valid).subtract(1))
                            .multiply(ee.Image(N_valid).multiply(2).add(5)).clip(roi);
  
  var n0 = ee.Image(N_valid).multiply(ee.Image(N_valid).subtract(1)).divide(2).clip(roi);
  
  var first_z = ee.List.repeat(zeros.clip(roi),4);
  
  var parameter_z = function(image, list_z){
    
    var n1 = ee.Image(ee.List(list_z).get(0));
    n1 = n1.add(image.multiply(image.subtract(1)).divide(2).clip(roi));
    list_z = ee.List(list_z).set(0,n1.rename('n1'));
    
    
    var m_xi = ee.Image(ee.List(list_z).get(1));
    m_xi = m_xi.add(image.multiply(image.subtract(1)).multiply(image.multiply(2).add(5)));
    list_z = ee.List(list_z).set(1,m_xi.rename('m_xi'));
    
    
    var m1 = ee.Image(ee.List(list_z).get(2));
    m1 = m1.add(image.multiply(image.subtract(1))
                 .divide(ee.Image(N).multiply(2).multiply(ee.Image(N).subtract(1))));
    list_z = ee.List(list_z).set(2,m1.rename('m1'));
    
    
    var m2 = ee.Image(ee.List(list_z).get(3));
    m2 = m2.add(image.multiply(image.subtract(1)).multiply(image.subtract(2))
    .divide(ee.Image(N).multiply(9).multiply(ee.Image(N).subtract(1)).multiply(ee.Image(N).subtract(2))));
    list_z = ee.List(list_z).set(3,m2.rename('m2'));
    
    return list_z;
  };
  
  var statistic_z = tie_xi.iterate(parameter_z,first_z);
  
  var m_xi = ee.List(statistic_z).get(1);
  
  var m1 = ee.List(statistic_z).get(2);
  
  var m2 = ee.List(statistic_z).get(3);
  
  var m = m0.subtract(m_xi).divide(18).add(m1).add(m2).sqrt();
  
  var cd = ee.Image(ee.List(tieCount).get(-3));
  var z = cd.divide(m);

  var trend_down = senslope.select(0).lt(0).and(z.abs().gt(1.960));
  var trend_up = senslope.select(0).gt(0).and(z.abs().gt(1.960));
  
  var trend = ones.where(trend_down,2).where(trend_up,3);
  
  return trend;
};

var PettittTest = function(dataSet,indicator,startYear){
  
  dataSet = dataSet.toList(dataSet.size());  
  
  var statistic_U = ee.List.sequence(0,dataSet.length().subtract(1)) 
    .map(function(t){
         var dataSet_xi = dataSet.slice(0,ee.Number(t).add(1));  
         var dataSet_xj = dataSet.slice(ee.Number(t).add(1),dataSet.length());
         var count = ee.Image.constant(0).clip(roi).set('ID','count');
    
         var accumulate_xi = function(image_xi,count){
             image_xi = image_xi.select(indicator);
             
             var accumulate_xj = function(image_xj,count){
                 image_xj = image_xj.select(indicator);
                 var diff = (image_xi.subtract(image_xj));
                 var sgn_diff = diff.where(diff.eq(0),0).where(diff.gt(0),1).where(diff.lt(0),-1);
                 count = ee.Image(count).add(sgn_diff.unmask());
                 return count;
             };
          
             count = ee.ImageCollection(dataSet_xj).iterate(accumulate_xj,count);
             
             return count;
         };
        
         count = ee.ImageCollection(dataSet_xi).iterate(accumulate_xi,count);
         //return ee.Image(count).abs();
         return ee.Image(count);//no absolute value here, but original value
    });    
  
 
  var K = ee.ImageCollection(statistic_U).max();//K+ used here 
  
  var P = K.pow(2).multiply(-6).divide(dataSet.length().pow(2).add(dataSet.length().pow(3)))
           .exp();//.multiply(2); //significance probability of K+
    
  var location = ee.Image.constant(-1000).clip(roi);  
  
  var index_K = function(image,location){
     var t = statistic_U.indexOf(image);
     var diff = ee.Image(image).subtract(K);
     location = ee.Image(location).where(diff.eq(0),t);
     return location;
  };
  
  var location_max = statistic_U.iterate(index_K,location);
  return ee.Image(location_max).add(startYear).addBands(P);
};  

//________________________________________________________________________________________
//                   data(******need to set time period manually******)
//________________________________________________________________________________________
var l8Col = dataPreprocessing("LANDSAT/LC08/C01/T1_SR","2000-1-1","2019-12-31");
var l7Col = dataPreprocessing("LANDSAT/LE07/C01/T1_SR","2000-1-1","2019-12-31");
var l5Col = dataPreprocessing("LANDSAT/LT05/C01/T1_SR","2000-1-1","2019-12-31");

//________________________________________________________________________________________
//                                 indicator calculation
//________________________________________________________________________________________
var l8Indicators = indicators(l8Col,'B2','B3','B4','B5','B6','B7');
//print(l8Indicators);
//Map.addLayer(l8Indicators);

var l7Indicators = indicators(l7Col,'B1','B2','B3','B4','B5','B7');
//print(l7Indicators);
//Map.addLayer(l7Indicators);

var l5Indicators = indicators(l5Col,'B1','B2','B3','B4','B5','B7');
//print(l5Indicators);
//Map.addLayer(l5Indicators);

var collection_merge = ee.ImageCollection(l8Indicators.merge(l7Indicators.merge(l5Indicators)));  
//print(collection_merge);
//Map.addLayer(collection_merge);

//________________________________________________________________________________________
//                          Annual mean indciator calculation 
//________________________________________________________________________________________
var annualMean_NDVI = annualMean(collection_merge,'NDVI',startYear,endYear,startMonth_NDVI,durationMonth_NDVI);
//print(annualMean_NDVI);
//Map.addLayer(annualMean_NDVI.select('NDVI_mean'));

var annualMean_AWEI_sh = annualMean(collection_merge,'AWEI_sh',startYear,endYear,1,12);
//print(annualMean_AWEI_sh);
//Map.addLayer(annualMean_AWEI_sh.select('AWEI_sh_mean'));

var annualMean_SWIR = annualMean(collection_merge,'SWIR',startYear,endYear,1,12);
//print(annualMean_SWIR);
//Map.addLayer(annualMean_SWIR.select('SWIR_mean'));

//________________________________________________________________________________________
//                          Time-series indicator preprocessing
//________________________________________________________________________________________
//join annualMean_AWEI_sh and annualMean_SWIR
var join_AM_AWEI_SWIR = ee.Join.inner().apply(annualMean_AWEI_sh,annualMean_SWIR, 
                            ee.Filter.equals({leftField: 'year',rightField: 'year'}));

var AWEI_sh_AMM = ee.ImageCollection(join_AM_AWEI_SWIR.map(function(feature) {
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'));
}));

//set awei<0 as -1, and awei>0 as 1, and awei as -1 when meanSWIR>0.1 and awei>0
var AWEI_sh_mean_modified = AWEI_sh_AMM.map(function(image) {
        
        var wl = image.select('AWEI_sh_mean');
        var sl = image.select('SWIR_mean');
        
        wl=wl.where(wl.lt(0),-1);
        wl=wl.where(wl.gt(0),1);
        
        //some urban pixels are misclassified as water, but their SWIR are higher than 0.1
        wl = wl.where(wl.eq(1).and(sl.gt(0.1)),-1);
        
        wl=image.addBands(wl.rename("AWEI_sh_mean_modified"));        
        wl=wl.select("AWEI_sh_mean_modified",'year');
        return wl;        
    });
    
//print(AWEI_sh_mean_modified);
//Map.addLayer(AWEI_sh_mean_modified.select('AWEI_sh_mean_modified'),{},'AWEI_sh');

//-----------------------------------------------------------------------------------------
//set NDVI<0 as 0, and NDVI>0.4 as 0.4
var NDVI_mean_modified =annualMean_NDVI.map(function(image) {
        var wl = image.select('NDVI_mean');
        wl=wl.where(wl.lt(0),0).where(wl.gt(0.4),0.4);        
        wl=image.addBands(wl.rename("NDVI_mean_modified"));        
        wl=wl.select("NDVI_mean_modified",'year');
        return wl;
    });
//print(NDVI_mean_modified);
//Map.addLayer(NDVI_mean_modified.select('NDVI_mean_modified'),{},'NDVI');

//________________________________________________________________________________________
//                              Long-term trend estimation
//________________________________________________________________________________________
var trend_NDVI = Trend(NDVI_mean_modified,'NDVI_mean_modified');

var trend_AWEI_sh = Trend(AWEI_sh_mean_modified,'AWEI_sh_mean_modified');

var trend_BCI = ee.Image("users/zhangyang2907/BCI_trend_FZ_0019");//(*****need to set manually*****)
//The website is: https://code.earthengine.google.com/?asset=users/zhangyang2907/BCI_trend_FZ_0019
       
var trend_all = trend_NDVI.multiply(100).add(trend_BCI.multiply(10))
                .add(trend_AWEI_sh).rename('classBand');
//print(trend_all);
//Map.addLayer(trend_all,{},'trend_all');

//_______________________________________________________________________________________________________
//      filtering urbanized pixels with abrupt changes from Class 111 
//_______________________________________________________________________________________________________

var PT_NDVI = PettittTest(NDVI_mean_modified,'NDVI_mean_modified',startYear).rename('cp_ndvi','P_ndvi');
//print(PT_NDVI);
//Map.addLayer(PT_NDVI,{},'PT_NDVI');

//extracting pixels with P<=0.5 and denotes them as 1112
var trend_1112_NDVI = trend_all.updateMask(trend_all.eq(111))
                         .where(PT_NDVI.select('P_ndvi').lte(0.5),1112);

//extracting the pixels with 1112 and the other pixels are masked                         
var trend_1112_NDVI_mask = trend_1112_NDVI.updateMask(trend_1112_NDVI.eq(1112));

var trend_all_new = trend_all.where(trend_1112_NDVI_mask.eq(1112),1112);
print(trend_all_new);
Map.addLayer(trend_all_new,{},'trend_all_new');

//________________________________________________________________________________________
//                              extracting urbanized classes
//________________________________________________________________________________________
var mask_urban = trend_all_new.eq(231).or(trend_all_new.eq(211)).or(trend_all_new.eq(131))
                          .or(trend_all_new.eq(332)).or(trend_all_new.eq(132))
                          .or(trend_all_new.eq(232)).or(trend_all_new.eq(212))  
                          .or(trend_all_new.eq(1112));
                          
var urban = trend_all_new.mask(mask_urban); 

//print(urban);
//Map.addLayer(urban,{},'urban');

//------------------------------------------------------------------------------------------------------
//                 Data preparation for eliminating falsely alarmed pixels
//------------------------------------------------------------------------------------------------------
var AMM_AWEI_List = AWEI_sh_mean_modified.toList(AWEI_sh_mean_modified.size());

var annualMean_SWIR_List = annualMean_SWIR.toList(annualMean_SWIR.size());

var annualMean_NDVI_List = annualMean_NDVI.toList(annualMean_NDVI.size());

//________________________________________________________________________________________
//                      filtering  water pixels 
//________________________________________________________________________________________

//using the mean AWEI of the last year to eliminate the water pixels
var meanAWEI_LY = ee.Image(AMM_AWEI_List.get(-1)).select('AWEI_sh_mean_modified');

//the pixels with meanAWEI_LY > 0 are denoted as 1
var urban_1 = urban.where(meanAWEI_LY.gt(0),1);
//print(urban_1);
//Map.addLayer(urban_1);

//________________________________________________________________________________________
//                       filtering wet pixels 
//________________________________________________________________________________________

var startYears = ee.List.sequence(0,19); //year t, t+1, and t+2 (*****need to set manually*****)
var min_3y_SWIR = ee.ImageCollection(startYears.map(function(t) {
        var min_3y = ee.ImageCollection(annualMean_SWIR_List.slice(t,ee.Number(t).add(3))).min(); 
        return min_3y;
    }));

//Using the max value of min_3y_SWIR.If the max value < 0.1,it is not possible to be urbanized pixels
var max_3yMin_SWIR = min_3y_SWIR.select('SWIR_mean').max();

//the pixels with max_AM_SWIR <0.1 are denoted as 2
var urban_2 = urban_1.where(max_3yMin_SWIR.lt(0.1).and(urban_1.neq(1)),2);
//print(urban_2);
//Map.addLayer(urban_2);

//_______________________________________________________________________________________________________
//                  filtering vegetated pixels 
//_______________________________________________________________________________________________________

var max_3y_NDVI = ee.ImageCollection(startYears.map(function(t) {
        var maxYear_3 = ee.ImageCollection(annualMean_NDVI_List.slice(t,ee.Number(t).add(3))).max(); 
        return maxYear_3;
    }));

//Using the min value of max_3y_NDVI.If the min value > 0.4,it is not possible to be urbanized pixels
var min_3yMax_NDVI = max_3y_NDVI.select('NDVI_mean').min();

//the pixels with min_3yMax_NDVI > 0.4 are denoted as 3
var urban_3 = urban_2
    .where(min_3yMax_NDVI.gt(0.4).and(urban_2.neq(1)).and(urban_2.neq(2)),3);
//print(urban_3);
//Map.addLayer(urban_3);

//_______________________________________________________________________________________________________
//                  filtering bare or urban pixels
//_______________________________________________________________________________________________________

var startYears_2 = ee.List.sequence(2,19); //year t, t-1, and t-2(*****need to set manually*****)

var min_3y_NDVI = ee.ImageCollection(startYears_2.map(function(t) {
        var minYear_3 = ee.ImageCollection
                        (annualMean_NDVI_List.slice(ee.Number(t).subtract(2),ee.Number(t).add(1))).min(); 
        return minYear_3;
    }));

//adding the first image (t=0) to min_3y_NDVI
min_3y_NDVI = ee.ImageCollection([annualMean_NDVI.first()]).merge(min_3y_NDVI);

//Using the max value of min_3y_NDVI.If the max value < 0.2,it is not possible to be urbanized pixels
var max_3yMin_NDVI = min_3y_NDVI.select('NDVI_mean').max();

//create a mask of water-related classes
var mask_waterUrban = trend_all_new.select('classBand').eq(332)
                               .or(trend_all_new.select('classBand').eq(132))
                               .or(trend_all_new.select('classBand').eq(232))
                               .or(trend_all_new.select('classBand').eq(212));
                               
//the pixels masked by max_3yMin_NDVI are denoted as 4
var urban_4 = urban_3.where(max_3yMin_NDVI.lt(0.2)
                     .and(urban_3.neq(1)).and(urban_3.neq(2)).and(urban_3.neq(3))
                     .and(mask_waterUrban.neq(1)),4);//only eliminate terrestrial pixels

//-------------------------------------------------------------------------------------------------
//       eliminating the vegetated pixels in cluster 1112
//-------------------------------------------------------------------------------------------------
var meanBCI_LY = ee.Image("users/zhangyang2907/meanBCI_LY_FZ_0019");//(*****need to set manually*****)
//The website is: https://code.earthengine.google.com/?asset=users/zhangyang2907/meanBCI_LY_FZ_0019

//the pixels with meanBCI_LY < 0 are denoted as 5 in the cluster 1112 
var urban_5 = urban_4.where(meanBCI_LY.lt(0)
                        .and(urban_4.neq(1)).and(urban_4.neq(2)).and(urban_4.neq(3))
                        .and(urban_4.neq(4)).and(urban_4.eq(1112)),5);
                        
//-------------------------------------------------------------------------------------------------
//       eliminating urban pixels with no change in cluster 1112
//-------------------------------------------------------------------------------------------------
// using mean SWIR of the first year to eliminate pixels in wet environment
var meanSWIR_FY = ee.Image(annualMean_SWIR_List.get(0)).select('SWIR_mean');

var min_3yMax_BCI = ee.Image("users/zhangyang2907/min_3yMax_BCI_FZ_0019");//(*****need to set manually*****)
//The website is: https://code.earthengine.google.com/?asset=users/zhangyang2907/min_3yMax_BCI_FZ_0019

//the pixels with min_3yMax_BCI > 0 are denoted as 6 in the cluster 1112 
var urban_6 = urban_5.where(min_3yMax_BCI.gt(0) //meaning BCI >0 in most time
                          .and(meanSWIR_FY.gt(0.1))  //meaning the land environment in the first year
                          .and(urban_5.neq(1)).and(urban_5.neq(2)).and(urban_5.neq(3))
                          .and(urban_5.neq(4)).and(urban_5.neq(5)).and(urban_5.eq(1112)),6);

print(urban_6);
Map.addLayer(urban_6,{},'urban_6');

//_______________________________________________________________________________________________________
//      ********************************   Output   **********************************************
//_______________________________________________________________________________________________________
//output of urban and stratified samples

//var sample = urban_6.stratifiedSample({numPoints:50,region:roi, scale:30,tileScale:16,geometries:true});
//-------------------------------------------------------------------------------------------------------
//get samples from a specific stratum
var sample = urban_6.updateMask(urban_6.eq(231))
                    .stratifiedSample({numPoints:85,region:roi, scale:30,tileScale:16,geometries:true});

//-------------------------------------------------------------------------------------------------------
//                                   Samples for water-related classes
//-------------------------------------------------------------------------------------------------------
//var urban_6_waterClusters = urban_6.eq(132).or(urban_6.eq(232)).or(urban_6.eq(332)).or(urban_6.eq(212));

//var sample = urban_6_waterClusters.updateMask(urban_6_waterClusters)
//                   .stratifiedSample({numPoints:30,region:roi, scale:30,tileScale:16,geometries:true});
//-------------------------------------------------------------------------------------------------------
print(sample);
Map.addLayer(sample,{},'sample');

Export.table.toDrive({
collection: sample,
description: 'sample',
});

Export.image.toDrive({
image: urban_6,
description: 'urban_6',
scale: 30,
region: roi
});

//________________________________________________________________________________________
//output of trend and stratified samples 

Export.image.toDrive({
image: trend_all_new,
description: 'trend_all_new',
scale: 30,
region: roi
});

//------------------------------------------------------------------------------------------------------
//                                  Samples for non-urbanized clusters
//------------------------------------------------------------------------------------------------------
//var sample_trend = trend_all_new.updateMask(trend_all_new.eq(121))
//                   .stratifiedSample({numPoints:60,region:roi, scale:30,tileScale:16,geometries:true});

//------------------------------------------------------------------------------------------------------
//                                     Samples for other class
//------------------------------------------------------------------------------------------------------
//var mask_other_urban = urban_6.eq(1).or(urban_6.eq(2)).or(urban_6.eq(3))
//                              .or(urban_6.eq(4)).or(urban_6.eq(5)).or(urban_6.eq(6));
                              
//var mask_other_trend = trend_all_new.neq(111).and(trend_all_new.neq(121))
//                        .and(trend_all_new.neq(311)).and(trend_all_new.neq(321))
//                        .and(trend_all_new.neq(231)).and(trend_all_new.neq(211))
//                        .and(trend_all_new.neq(131)).and(trend_all_new.neq(332))
//                        .and(trend_all_new.neq(132)).and(trend_all_new.neq(232))
//                        .and(trend_all_new.neq(212)).and(trend_all_new.neq(1112));
                        
//var mask_other = trend_all_new.updateMask(mask_other_trend).unmask(mask_other_urban).neq(0);

//var sample_trend = mask_other.updateMask(mask_other)
//                   .stratifiedSample({numPoints:85,region:roi, scale:30,tileScale:16,geometries:true});
//--------------------------------------------------------------------------------------------------------
//print(sample_trend);
//Map.addLayer(sample_trend,{},'sample_trend');

//Export.table.toDrive({
//collection: sample_trend,
//description: 'sample_trend',
//});






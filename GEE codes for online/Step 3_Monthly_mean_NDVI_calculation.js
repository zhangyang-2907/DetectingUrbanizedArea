//----------------------------------------------------------------------------------------
//                      Step 3: calculating monthly mean NDVI
//----------------------------------------------------------------------------------------
//Note: there are one place that need to set the temporal period manually. 
//It is in Data section (Lines 77-79).
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

//(2) Data Preprocessing
var dataPreprocessing = function (dataName,startDate,endDate){
    var landsat_collection = ee.ImageCollection(dataName);
    var landsat_studyArea = landsat_collection.filterBounds(roi);
    var landsat = landsat_studyArea.filterDate(startDate,endDate);
    var landsat_imageNoClouds = landsat.map(maskClouds); 
    return landsat_imageNoClouds.map(function(img){return img.clip(roi);});
};

//(3) NDVI calculatioin
var indicators = function (dataSet,BLUE,GREEN,RED,NIR,SWIR1,SWIR2){   
    dataSet = dataSet.map(function(image) {
             var ndvi = image.normalizedDifference([NIR, RED]);
             ndvi = ndvi.where(ndvi.lt(0),0);
             return  image.addBands(ndvi.rename("NDVI")).toFloat();
          });
          
    return dataSet;
};

//(4) Monthly mean NDVI calculatioin
var monthlyMean = function (dataSet, indicator){
  
    var months = ee.List.sequence(1, 12);
    var reducers = ee.Reducer.mean().combine({reducer2: ee.Reducer.stdDev(),sharedInputs: true});
    
    var monthlyMean = ee.ImageCollection.fromImages(
      months.map(function (m) {
        var monthlyMean_indi = dataSet.select(indicator).filter(ee.Filter.calendarRange(m, m, 'month'))
                                      .mean()
                                      .set('month', m);
                                      
        // calculate the mean NDVI values of each month for the study area                              
        var monthlyMean_stdDev = monthlyMean_indi.reduceRegion(
                           {reducer:reducers, geometry: roi, scale:30, maxPixels: 1e13, tileScale:16});
        
        return monthlyMean_indi
               .set('mean', monthlyMean_stdDev.get('NDVI_mean'),'SD', monthlyMean_stdDev.get('NDVI_stdDev'));
      }));
   
    return monthlyMean;
};


//________________________________________________________________________________________
//                                      Implementation
//________________________________________________________________________________________

//(1) Data (******need to set the temporal periods manually******)
var l8Col = dataPreprocessing("LANDSAT/LC08/C01/T1_SR","2000-1-1","2019-12-31");
var l7Col = dataPreprocessing("LANDSAT/LE07/C01/T1_SR","2000-1-1","2019-12-31");
var l5Col = dataPreprocessing("LANDSAT/LT05/C01/T1_SR","2000-1-1","2019-12-31");
//print(l8Col);
//print(l7Col);
//print(l5Col);

//(2) Calculating NDVI
var l8Indicators = indicators(l8Col,'B2','B3','B4','B5','B6','B7');
var l7Indicators = indicators(l7Col,'B1','B2','B3','B4','B5','B7');
var l5Indicators = indicators(l5Col,'B1','B2','B3','B4','B5','B7');
//print(l8Indicators);
//print(l7Indicators);
//print(l5Indicators);
var collection_merge = ee.ImageCollection(l8Indicators.merge(l7Indicators.merge(l5Indicators)));  
//print (collection_merge);

//(3) Calculating monthly mean NDVI
var monthlyMean_NDVI = monthlyMean(collection_merge,'NDVI');
print(monthlyMean_NDVI);
Map.addLayer(monthlyMean_NDVI);

//(4) create chart for showing monthly mean NDVI in which the stable period can be observed and determined
var trend_mean = monthlyMean_NDVI.toList(monthlyMean_NDVI.size());
  var array_x = ee.List([]);
  var array_y = ee.List([]);
  var array_SD = ee.List([]);
  var first = ee.List([array_x, array_y, array_SD]);
  
  var Array_indi_mean = function(image, list){
     var added_mean = ee.Image(image).get('mean');
     var added_stdDev = ee.Image(image).get('SD');
     var added_year = ee.Image(image).get('month');
     
     var array_x = ee.List(ee.List(list).get(0)).add(added_year);
     var array_y = ee.List(ee.List(list).get(1)).add(added_mean);
     var array_SD = ee.List(ee.List(list).get(2)).add(added_stdDev);
     list = ee.List([array_x, array_y, array_SD]);
     return list;
  };
  
  var Array_forChart = trend_mean.iterate(Array_indi_mean, first);
  array_x = ee.List(Array_forChart).get(0);
  array_y = ee.List(Array_forChart).get(1);
  array_SD = ee.List(Array_forChart).get(2);
  
  var yValues = ee.Array.cat([array_y, array_SD], 1);
  
  var chart = ui.Chart.array.values(yValues, 0, array_x);
  print(chart);
  
 











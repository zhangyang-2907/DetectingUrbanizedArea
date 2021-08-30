//----------------------------------------------------------------------------------------
//                      Step 2: calculating BCI trend
//----------------------------------------------------------------------------------------
//Note: there are some places that need to set data sources manually (Lines 33 and 37)

//var BCI_1 = ee.Image("users/zhangyang2907/BCI_FZ_0009");
//The website is: https://code.earthengine.google.com/?asset=users/zhangyang2907/BCI_FZ_0009

//var BCI_2 = ee.Image("users/zhangyang2907/BCI_FZ_1019");
//The website is: https://code.earthengine.google.com/?asset=users/zhangyang2907/BCI_FZ_1019

//----------------------------------------------------------------------------------------
//Note: there are some places that need to set list sequence values manually(Lines 49 and 72)

//var BCI_collection = ee.ImageCollection(ee.List.sequence(0,39,2).map(function(t){

//var startYears_2 = ee.List.sequence(2,19); 

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
//             Merging BCI images of different periods into one image
//________________________________________________________________________________________
var BCI_1 = ee.Image("users/zhangyang2907/BCI_FZ_0009");
BCI_1 = BCI_1.select(BCI_1.bandNames().remove('constant'));
//print(BCI_1);

var BCI_2 = ee.Image("users/zhangyang2907/BCI_FZ_1019");
BCI_2 = BCI_2.select(BCI_2.bandNames().remove('constant'));
//print(BCI_2);

var BCI_3 = BCI_1.addBands(BCI_2);
//print(BCI_3);
//Map.addLayer(BCI_3);

//________________________________________________________________________________________
//              Converting BCI image bands into an imagecollection 
//            (******need to set list sequence values manually******)
//________________________________________________________________________________________
var BCI_collection = ee.ImageCollection(ee.List.sequence(0,39,2).map(function(t){
  var image = BCI_3.slice(t,ee.Number(t).add(2)).rename('BCI_mean','year');
  return image;
}));

print(BCI_collection);
Map.addLayer(BCI_collection.select('BCI_mean')); 

//________________________________________________________________________________________
//                        Calculating meanBCI_lastYear
//________________________________________________________________________________________

//using the mean BCI of the last year to eliminate the vegetation pixels
var AMM_BCI_List = BCI_collection.toList(BCI_collection.size());
var meanBCI_lastYear = ee.Image(AMM_BCI_List.get(-1)).select('BCI_mean').rename('MeanBCI_lastYear');

//print(meanBCI_lastYear);
//Map.addLayer(meanBCI_lastYear,{},'meanBCI_lastYear');

//________________________________________________________________________________________
//                       Calculating 3-year maximum BCI
//          (******need to set list sequence values manually******)
//________________________________________________________________________________________
var startYears_2 = ee.List.sequence(2,19); //year t, t-1, and t-2

// 3 year maximum of BCI between t and t-2
var maxBCI_3year = ee.ImageCollection(startYears_2.map(function(t) {
        var maxYear_3 = ee.ImageCollection
                        (AMM_BCI_List.slice(ee.Number(t).subtract(2),ee.Number(t).add(1))).max(); 
        return maxYear_3;
    }));

//adding the first image (t=0) to maxBCI_3year
maxBCI_3year = ee.ImageCollection([BCI_collection.first()]).merge(maxBCI_3year);
//print(maxBCI_3year);
//Map.addLayer(maxBCI_3year.select('BCI_mean'));

//Using the min value of maxBCI_3year.If the min value > 0,it is not possible to be urbanized pixels
//because there is no vegetation converted into bare land or impervious surfaces
var min_3yMax_BCI = maxBCI_3year.select('BCI_mean').min().rename('min_3yMax_BCI');
//print(min_3yMax_BCI);
//Map.addLayer(min_3yMax_BCI,{},'min_3yMax_BCI');

//________________________________________________________________________________________
//                       calculating the trend of BCI
//________________________________________________________________________________________

//set BCI<0 as 0
var BCI_mean_modified =BCI_collection.map(function(image) {
        var wl = image.select('BCI_mean');
        wl=wl.where(wl.lt(0),0);        
        wl=image.addBands(wl.rename("BCI_mean_modified"));        
        wl=wl.select("BCI_mean_modified",'year');
        return wl;
    });


var Trend = function(dataSet,indicator){
  
  var senslope = dataSet.reduce(ee.Reducer.sensSlope());
 
  var N_valid = dataSet.select(indicator).count();
 
  var ones = ee.Image.constant(1).clip(roi);   
  var zeros = ee.Image.constant(0).clip(roi); 
  dataSet = dataSet.toList(dataSet.size());
  var N = dataSet.length();
  
  var countArray = ee.List.sequence(0,N.subtract(1))
    .map(function (xi) {return ee.Image.constant(-1).clip(roi).set('a7',xi);});
  
  countArray = countArray.add(zeros);
  countArray = countArray.add(zeros);
  countArray = countArray.add(zeros);
  
 
  var first = countArray.add(ones);
  first = first.add(zeros);
 
  var dataSet_xi = dataSet.slice(0,N.subtract(1));
 
  var a1 = function(image,list){
          
          var xi = ee.List(dataSet).indexOf(image);
          var a2 = image.select(indicator);
          var first_xj = ee.List(list).set(-1, a2);
          
          var a3 = dataSet.slice(ee.Number(xi).add(1),dataSet.length());
          
          var a4 = function(a5,a6){
             var xj = ee.List(dataSet).indexOf(a5);
             var currentImage = a5.select(indicator);
             var previousImage = ee.Image(ee.List(a6).get(-1));
             var diff = (currentImage.subtract(previousImage));
            
             var a7 = ee.Image(ee.List(a6).get(xi));
             var a8 = ee.Image(ee.List(a6).get(xj));
             var count = ee.Image(ee.List(a6).get(-2));
             
             var a9 = diff.eq(0).and(a7.neq(0));
             a8 = a8.where(a9, zeros);
             count = count.where(a9,count.add(1));
             a6 = ee.List(a6).set(xj,a8);
             a6 = ee.List(a6).set(-2,count);
            
             var c_diff = diff.where(diff.gt(0),1).where(diff.neq(1),0);
             var a10 = ee.Image(ee.List(a6).get(-5));
             a10 = a10.add(c_diff.unmask());
             a6 = ee.List(a6).set(-5,a10);
             
             var d_diff = diff.where(diff.eq(0),0).where(diff.gt(0),0).where(diff.lt(0),1);
             var a11 = ee.Image(ee.List(a6).get(-4));
             a11 = a11.add(d_diff.unmask());
             a6 = ee.List(a6).set(-4,a11);
             
             var cd_diff = diff.where(diff.eq(0),0).where(diff.gt(0),1).where(diff.lt(0),-1);
             var cd = ee.Image(ee.List(a6).get(-3));
             cd = cd.add(cd_diff.unmask());
             a6 = ee.List(a6).set(-3,cd);
    
             return ee.List(a6);
             
          };
          list = ee.ImageCollection(a3).iterate(a4,first_xj);
          
          list = ee.List(list).set(xi,ee.List(list).get(-2));
          list = ee.List(list).set(-2,ee.Image.constant(1).clip(roi));
          return ee.List(list);
  };
  var tieCount = ee.ImageCollection(dataSet_xi).iterate(a1,first);
  
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

var trend_BCI = Trend(BCI_mean_modified,'BCI_mean_modified');

print(trend_BCI);
Map.addLayer(trend_BCI);

//________________________________________________________________________________________
//********************************   Output   ********************************************
//________________________________________________________________________________________

Export.image.toAsset({
  image: trend_BCI,
  description: 'BCI_trend_FZ_0019',
  scale: 30,  
  region: roi
});
  
Export.image.toAsset({
  image: meanBCI_lastYear,
  description: 'MeanBCI_LY_FZ_0019',
  scale: 30,  
  region: roi
});

Export.image.toAsset({
  image: min_3yMax_BCI,
  description: 'min_3yMax_BCI_FZ_0019',
  scale: 30,  
  region: roi
});









  
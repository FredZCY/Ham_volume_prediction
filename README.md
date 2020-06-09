# Ellipsoidal Ham Volume Prediction

This work determine the volume of hams particularly for the horizontal viewpoint. Three mathematical approaches had been proposed to detect the defective images in the image sequences:
1) Based on the k-nearest neighbor
2) Based on minor axis
3) Based on Y-direction

The method had been evaluated on 2 databases: HamDB-A and HamDB-B, where the data elicitation method are slightly different.
Note that, a Mask Region–based convolutional neural network (Mask R-CNN) object segmentation approach was adopted to extract the volume of the ham object from a video.

Software is written and tested using Matlab 2020a, toolbox required:

Statistics and Machine Learning Toolbox

The files include:
1) defectDetection_knn.m - defect detection based on the knn 
2) defectDetection_minor.m - defect detection based on the minor axis
3) defectDetection_y_axis.m - defect detection based on the y direction
4) curveFitting_polynomial.m - optimal curve fitting functions based on polynomial order 2
5) curveFitting_power.m - optimal curve fitting functions based on power function



Acquisition setup for the elicitation and video recording of the ham sample:
<img src="https://github.com/christy1206/Ham_volume_prediction/blob/pic/data_elicitation.JPG" width="500" height="300"/>



Example of 16 ham samples in HAMDB-A:

<img src="https://github.com/christy1206/Ham_volume_prediction/blob/pic/HamDB-A.JPG" width="700" height="500"/>



Example of 16 ham samples in HAMDB-B:

<img src="https://github.com/christy1206/Ham_volume_prediction/blob/pic/HamDB-B.JPG" width="700" height="500"/>



(a) Good Mask R-CNN image and (b) Poor Mask R-CNN image:

<img src="https://github.com/christy1206/Ham_volume_prediction/blob/pic/maskrcnn.JPG" width="500" height="200"/>


Volume prediction result on HamDB-A (left) and HamDB-B (right):

<img src="https://github.com/christy1206/Ham_volume_prediction/blob/pic/result1.JPG" width="300" height="300"/> <img src="https://github.com/christy1206/Ham_volume_prediction/blob/pic/result2.JPG" width="300" height="300"/>




If you use this method in your research, please cite:

@article{liong2020ham,\
title={A Statistical Approach in Enhancing the Volume Prediction of Ellipsoidal Ham},\
author={Y. S. Gan and Lan Wei and Yiming Hanb and Chenyu Zhang and Yen-Chang Huang and Sze-Teng Liong},\
journal={Journal of Food Engineering},\
volume={1},\
pages={1--1},\
year={2020},\
publisher={Elsevier}\
}

If you have suggestions or questions regarding this method, please reach out to stliong@fcu.edu.tw

Thank you for your interest and support.

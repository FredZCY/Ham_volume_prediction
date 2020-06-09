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




Acquisition setup for the elicitation and video recording of the ham sample.

Example of 16 ham samples in HAMDB-A

Example of 16 ham samples in HAMDB-B

(a) Good Mask R-CNN image and (b) Poor Mask R-CNN image

If you use this method in your research, please cite:

@article{liong2020ham,
title={A Statistical Approach in Enhancing the Volume Prediction of Ellipsoidal Ham},
author={Y. S. Gan and Lan Wei and Yiming Hanb and Chenyu Zhang and Yen-Chang Huang and Sze-Teng Liong},
journal={Journal of Food Engineering},
volume={1},
pages={1--1},
year={2020},
publisher={Elsevier}
}

If you have suggestions or questions regarding this method, please reach out to stliong@fcu.edu.tw

Thank you for your interest and support.

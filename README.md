# Investigating Trajectories of Actin-based Motility

This is the independent work done by Zach Lyu as a research assistant at Dr. Yuan Lin's Lab.

After the Listeria “highjacked” host
cells, It was then surrounded by actin
filaments. As it grew, actin filaments
began to form a “comet tail”.

Despite previous research, questions remain:
(1)	What are the factors regulating the elongating direction and force generation of the actin network? 
(2)	How does the load geometry/size affect actin-based motion? 
(3)	What is the full motion of actin-propelled load? What are the implications for force generation mechanism?
(4)	What causes the dramatic morphological transformation of the nuclear envelope during closed mitosis of fission yeast Schizosaccharomyces pombe?

In this project, two algorithms were designed to analyze actin-based mitility trajectories in 2D and 3D respectively.

## 2D

1. Open a png file (BW curve after processed by imageJ).
![sample 1](https://github.com/bijiuni/actin_motility/blob/master/2D/1.png)
![sample 2](https://github.com/bijiuni/actin_motility/blob/master/2D/2.png)

2. Choose a first point as origin (not all the parts of the curve should be used to calculate curvature, a part of the curve is chosen) and the whole curve becomes green if successful. (Otherwise the chosen point is not on the curve and the the whole code ends)
![sample 1](https://github.com/bijiuni/actin_motility/blob/master/2D/1green.jpg)
![sample 2](https://github.com/bijiuni/actin_motility/blob/master/2D/2green.jpg)

3. Choose the second point as ending point (and the chosen curve will become red, otherwise screen display 'not on curve')
![sample 1](https://github.com/bijiuni/actin_motility/blob/master/2D/1red.jpg)
![sample 2](https://github.com/bijiuni/actin_motility/blob/master/2D/2red.jpg)

4. If all of above is successful, the code randomly picks 10 points from the chosen part of the curve and calculate curvature by fitting polygons to points and display the results.


## 3D

The case in 3D is more complicated.

First the 3D model needs to be smoothed. The resulting samples are shown below.
[Sample 1](https://github.com/bijiuni/actin_motility/blob/master/3D/Smoothed-2.tif)
[Sample 2](https://github.com/bijiuni/actin_motility/blob/master/3D/Smoothed3.tif)


Then the curvatures need to be skeletonized. This is done by MATLAB code and Image J.
![sample 1](https://github.com/bijiuni/actin_motility/blob/master/3D/Figure%201.jpg)
![sample 2](https://github.com/bijiuni/actin_motility/blob/master/3D/Figure%202.jpg)


Finally the codes calculate the curvature and output the result:
![sample output](https://github.com/bijiuni/actin_motility/blob/master/3D/Curvature%20answers.jpg)


## Built With

* [ImageJ](https://www.tensorflow.org/) - An open source machine learning framework
* [Matlab](https://www.mathworks.com/products/matlab.html)

## Authors

* **Zach Lyu** - [bijiuni](https://github.com/bijiuni)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


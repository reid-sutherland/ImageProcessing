#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

//Opencv
#include "opencv2/opencv.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace std;
using namespace cv;

// **************************************************************
// Function for finding MSE between two matrices
// **************************************************************
float findMSE(Mat a, Mat b) {
	// assert matching sizes
	if (a.size != a.size) {
		cout << "Matrices are not the same size, cannot compute MSE, beep boop" << endl;
		return -1;
	}

	// *********************************************
	// MSE = (1/N) * \sigma(val - val')^2
	// MSE = mean of all squared errors
	// *********************************************

	float N = (float) a.rows * a.cols;
	float sigmaN = 0.0f;
	float MSE = 0.0f;

	for (int j = 0; j < a.rows; j++) {
		for (int i = 0; i < a.cols; i++) {

			// error
			float error = (float)a.at<uchar>(i, j) - (float)b.at<uchar>(i, j);

			// squared error
			error = powf(error, 2.0f);

			// add squared error to sum
			sigmaN += error;
		}
	}

	// get mean from sum
	MSE = sigmaN / N;

	// print
	cout << "MSE: " << MSE << endl;
	// cout << "sigmaN: " << sigmaN << endl;
	cout << endl;

	return MSE;
}


// **************************************************************
// Laplacian Filters
// **************************************************************
void laplacian() {

	Mat original = imread("../images/megan.jpg", CV_32F);

	Mat img;

	cvtColor(original, img, COLOR_BGR2GRAY);

	float values1[9] = {1, 0, 1, 0, -4, 0, 1, 0, 1};
	Mat filter1(3, 3, CV_32F, values1);

	float values2[9] = {1, 1, 1, 1, -8, 1, 1, 1, 1};
	Mat filter2(3, 3, CV_32F, values2);

	float values3[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
	Mat sobel1(3, 3, CV_32F, values3);

	float values4[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
	Mat sobel2(3, 3, CV_32F, values4);

	float valuesD[25] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -24, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	Mat filterD(5, 5, CV_32F, valuesD);


	Mat img1(img.rows, img.cols, CV_32F);
	Mat img2(img.rows, img.cols, CV_32F);
	Mat imgSobel1(img.rows, img.cols, CV_32F);
	Mat imgSobel2(img.rows, img.cols, CV_32F);
	Mat imgD(img.rows, img.cols, CV_32F);

	filter2D(img, img1, -1, filter1);
	filter2D(img, img2, -1, filter2);
	filter2D(img, imgSobel1, -1, sobel1);
	filter2D(img, imgSobel2, -1, sobel2);
	filter2D(img, imgD, -1, filterD);

	Mat img2T = img2.clone();
	for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.cols; j++) {
			if (img2T.at<uchar>(i, j) < 100)
				img2T.at<uchar>(i, j) = 0;
			else
				img2T.at<uchar>(i, j) = 255;
		}
	}


	imshow("original", original);
	imshow("gray", img);
	imshow("filter 1", img1);
	imshow("filter 2", img2);
	imshow("filter 2 T", img2T);
	// imshow("sobel 1", imgSobel1);
	// imshow("sobel 2", imgSobel2);
	// imshow("filter D", imgD);
	while(waitKey(30) != 27);
	destroyAllWindows();

}


// **************************************************************
// Laplacian Filters
// **************************************************************
void highPassEqualization() {

	Mat img = imread("lena.tif", CV_8U);
	// Mat img = imread("Bean.jpg", CV_8U);

	Mat imgFilter(img.rows, img.cols, CV_8U);
	Mat imgFilterEq(img.rows, img.cols, CV_8U);
	Mat imgEq(img.rows, img.cols, CV_8U);
	Mat imgEqFilter(img.rows, img.cols, CV_8U);
	Mat imgAvg(img.rows, img.cols, CV_8U);
	Mat imgLap(img.rows, img.cols, CV_8U);

	cout << "Sizeof image: " << img.size() << endl;
	int numPixels = img.rows * img.cols;

	//! High Pass Filter
	float a = 0;
	float b = -1.0;
	float c = 9.0;
	float values[9] = {a, b, a, b, c, b, a, b, a};
	//
	// Mat highPassFilter(3, 3, CV_32F, values);
	// // scale filter by 1/4
	// highPassFilter *= (1.0  / 4);
	// cout << highPassFilter << endl;
	//
	// // filter then equalize first
	// filter2D(img, imgFilter, -1, highPassFilter);
	// equalizeHist(imgFilter, imgFilterEq);
	//
	// // now equalize then filter
	// equalizeHist(img, imgEq);
	// filter2D(imgEq, imgEqFilter, -1, highPassFilter);


	//! Average Filter
	a = 1.0;
	b = 0.0;
	c = 1.0;
	float values2[9] = {a, b, a, b, c, b, a, b, a};

	Mat avgFilter(3, 3, CV_32F, values2);
	// scale filter by 1/5
	avgFilter *= (1.0 / 5);
	cout << "avgFilter\n";
	cout << avgFilter << endl << endl;

	// filter
	filter2D(img, imgAvg, -1, avgFilter);


	//! Laplacian Filter
	a = 1.0;
	b = 0.0;
	c = -4.0;
	float values3[9] = {a, b, a, b, c, b, a, b, a};

	Mat lapFilter(3, 3, CV_32F, values3);
	// // scale filter by 1/4
	// lapFilter *= (1.0 / 4);
	cout << "lapFilter\n";
	cout << lapFilter << endl << endl;

	// filter
	filter2D(img, imgLap, -1, lapFilter);



	// show images
	if (!img.empty())
	imshow("original", img);
	if (!imgAvg.empty())
	imshow("average", imgAvg);
	if (!imgLap.empty())
	imshow("laplacian", imgLap);
	// if (!imgFilter.empty())
	// imshow("filter", imgFilter);
	// if (!imgFilterEq.empty())
	// imshow("filterEq", imgFilterEq);
	// if (!imgEq.empty())
	// imshow("eq", imgEq);
	// if (!imgEqFilter.empty())
	// imshow("eqFilter", imgEqFilter);
	while(waitKey(30) != 27);
	destroyAllWindows();

	// write images
	imwrite("../imagesMidterm/img.png", img);
	imwrite("../imagesMidterm/imgAvg.png", imgAvg);
	imwrite("../imagesMidterm/imgLap.png", imgLap);
	// imwrite("../imagesMidterm/filter.png", imgFilter);
	// imwrite("../imagesMidterm/eq.png", imgEq);
	// imwrite("../imagesMidterm/filterEq.png", imgFilterEq);
	// imwrite("../imagesMidterm/eqFilter.png", imgEqFilter);
}

// **************************************************************
// Problem Six - Histograms and Equalization
// **************************************************************
int problemSix() {
	// bean image
	Mat bean = imread("Bean.jpg", CV_8U);
	// create data ptr for bean
	// uchar* beanData = bean.data;

	// gray histogram - 256 gray valueslaplaceC
	Mat hist = Mat::zeros(1, 256, CV_32F);

	// width should be twice the number of values
	int histWidth = 2 * hist.cols;
	int histHeight = 2 * hist.cols;

	// mat for displaying histogram, set height to max gray count
	Mat histImg = Mat::zeros(histHeight, histWidth, CV_8U);

	// create histogram
	for (int y = 0; y < bean.rows; y++) {
		for (int x = 0; x < bean.cols; x++) {
			uchar grayVal = bean.data[y*bean.cols + x];
			// increment the histogram col for the gray value
			hist.at<float>((int)grayVal)++;
		}
	}

	// find max value of hist
	double maxVal;
	minMaxLoc(hist, NULL, &maxVal, NULL, NULL);

	// draw the hist chart
	for (int i = 1; i < hist.cols; i++) {

		int scaledStartVal = cvRound((int) ( (hist.at<float>(i-1) / maxVal) * histHeight ));
		int scaledEndVal = cvRound((int) ( (hist.at<float>(i) / maxVal) * histHeight ));

		int pStartX = 2 * (i - 1);
		int pStartY = histHeight - scaledStartVal;
		int pEndX = 2 * i;
		int pEndY = histHeight - scaledEndVal;

		line(histImg, Point(pStartX, pStartY), Point(pEndX, pEndY), Scalar(255, 255, 255));
	}

	imshow("Our Histogram", histImg);
	while(waitKey(30) != 27);
	destroyAllWindows();


	// now calculate histogram using opencv's calcHist
	Mat opencvHist;
	int channels[] = {0};
	int histSize = hist.cols;
	float range[] = {0, 256};
	const float* histRange[] = {range};

	calcHist(&bean, 1, channels, Mat(), opencvHist, 1, &histSize, histRange);

	double maxValue;
	minMaxLoc(opencvHist, NULL, &maxValue, NULL, NULL);

	Mat opencvHistImg = Mat::zeros(histHeight, histWidth, CV_8U);

	// draw the histogram
	for (int i = 1; i < histSize; i++) {
		int scaledStartVal = cvRound((int) ( (opencvHist.at<float>(i-1) / maxValue) * histHeight ));
		int scaledEndVal = cvRound((int) ( (opencvHist.at<float>(i) / maxValue) * histHeight ));

		int pStartX = 2 * (i - 1);
		int pStartY = histHeight - scaledStartVal;
		int pEndX = 2 * i;
		int pEndY = histHeight - scaledEndVal;

		line(opencvHistImg, Point(pStartX, pStartY), Point(pEndX, pEndY), Scalar(255, 255, 255));
	}

	// show
	imshow("Opencv Histogram", opencvHistImg);
	while(waitKey(30) != 27);
	destroyAllWindows();

	// write
	imwrite("../images(#6)/histImg.png", histImg);
	imwrite("../images(#6)/opencvHistImg.png", opencvHistImg);


	// *********************************************
	// Equalize the histogram
	// *********************************************

	// cumulative density function
	Mat cdf = Mat::zeros(1, 256, CV_32F);
	cdf.at<float>(0) = hist.at<float>(0);
	for (int i = 1; i < hist.cols; i++) {
		float val = hist.at<float>(i);
		cdf.at<float>(i) = val + cdf.at<float>(i-1);
	}

	// cumulative probability function
	float numPixels = cdf.at<float>(cdf.cols-1);
	Mat cpf = Mat::zeros(1, 256, CV_32F);
	for (int i = 0; i < hist.cols; i++) {
		cpf.at<float>(i) = cdf.at<float>(i) / numPixels;
	}

	// scale cpf to new image
	Mat histEq(cpf.rows, cpf.cols, CV_8U);
	for (int i = 0; i < hist.cols; i++) {
		histEq.at<uchar>(i) = (uchar) floorf(cpf.at<float>(i) * 255);
	}

	Mat beanEq(bean.rows, bean.cols, CV_8U);
	for (int i = 0; i < bean.cols; i++) {
		for (int j = 0; j < bean.rows; j++) {
			beanEq.at<uchar>(i, j) = histEq.at<uchar>((int) bean.at<uchar>(i, j));
		}
	}

	// draw the histogram Equalized
	Mat histEqImg(bean.rows, bean.cols, CV_8U);
	for (int i = 1; i < histEq.cols; i++) {
		int scaledStartVal = cvRound((int) ( (histEq.at<float>(i-1) / maxValue) * histEq.rows ));
		int scaledEndVal = cvRound((int) ( (histEq.at<float>(i) / maxValue) * histEq.rows ));

		int pStartX = 2 * (i - 1);
		int pStartY = histEq.rows - scaledStartVal;
		int pEndX = 2 * i;
		int pEndY = histEq.rows - scaledEndVal;

		line(histEqImg, Point(pStartX, pStartY), Point(pEndX, pEndY), Scalar(255, 255, 255));
	}


	// *********************************************
	// Compare to OpenCV's EqualizeHist
	// *********************************************

	Mat opencvBeanEq;
	equalizeHist(bean, opencvBeanEq);
	Mat opencvBeanEq2;
	equalizeHist(opencvBeanEq, opencvBeanEq2);

	// write
	imwrite("../images(#6)/histEq.png", histEq);
	imwrite("../images(#6)/histEqImg.png", histEqImg);
	imwrite("../images(#6)/opencvBeanEq.png", opencvBeanEq);
	imwrite("../images(#6)/opencvBeanEq2.png", opencvBeanEq2);


	// show
	imshow("Equalized Bean Image", opencvBeanEq);
	imshow("Equalized (x2) Bean Image", opencvBeanEq2);
	// imshow("Equalized Bean Image", beanEq);
	imshow("Original Bean Image", bean);
	while(waitKey(30) != 27);
	destroyAllWindows();


	return 0;
}



// **************************************************************
// Problem Seven - Image Quantization and Smoothing
// **************************************************************

int problemSeven() {

	Mat lena = imread("lena.tif", CV_8U);
	imwrite("../images(#7)/original.png", lena);

	Mat lenaNoisy;
	Mat quan(lena.rows, lena.cols, CV_8U);
	Mat noisy(lena.rows, lena.cols, CV_8U);
	Mat noisyQuan(lena.rows, lena.cols, CV_8U);
	Mat quanFiltered(lena.rows, lena.cols, CV_8U);
	Mat noisyFiltered(lena.rows, lena.cols, CV_8U);
	Mat noisyQuanFiltered(lena.rows, lena.cols, CV_8U);

	// *********************************************
	// Quantize the image
	// *********************************************
	for (int j = 0; j < lena.rows; j++) {
		for (int i = 0; i < lena.cols; i++) {

			// get gray value and cast to int
			int grayVal = (int) lena.at<uchar>(i, j);

			// 256 divided by 8 is 32, so we need to take each value,
			// compute modulo 32, then compare that to the middle of each
			// layer (values where modulo 32 == 16), and adjust the
			// value accordingly
			int mod32 = grayVal % 32;
			int difference = 16 - mod32;
			int newVal = grayVal + difference;

			// assign quantized value to quantized image
			quan.at<uchar>(i, j) = (uchar) newVal;

			// print
			// cout << "lena(" << i << ", " << j << "): " << (int) lena.at<uchar>(i, j);
			// cout << "\tlenaQuan(" << i << "," << j << "): " << (int) quan.at<uchar>(i, j) << endl;
			// cout << endl;
		}
	}

	cout << "*** Mean Squared Error between image and quantized image ***" << endl;
	findMSE(lena, quan);

	// show
	imshow("Original", lena);
	imshow("Quantized", quan);
	while(waitKey(30) != 27);
	destroyAllWindows();

	// write
	imwrite("../images(#7)/quantized.png", quan);



	// *********************************************
	// Generate random noise
	// *********************************************

	// random seed
	srand(time(NULL));

	for (int j = 0; j < lena.rows; j++) {
		for (int i = 0; i < lena.cols; i++) {

			// get random noise - from -16 to 16
			int noise = rand() % 33 - 16;

			// add noise to lena
			noisy.at<uchar>(i, j) = lena.at<uchar>(i, j) + (uchar) noise;

			// print
			// if (j < 2 && i < 50) {
			// 	cout << "lena(" << i << ", " << j << "): " << (int) lena.at<uchar>(i, j);
			// 	cout << "\tlenaNoisy(" << i << ", " << j << "): " << (int) noisy.at<uchar>(i, j) << endl;
			// 	cout << endl;
			// }
		}
	}

	cout << "*** Mean Squared error between image and noisy image ***" << endl;
	findMSE(lena, noisy);

	// show
	imshow("Original", lena);
	imshow("Quantized", quan);
	imshow("Noisy", noisy);
	while(waitKey(30) != 27);
	destroyAllWindows();

	// write
	imwrite("../images(#7)/noisy.png", noisy);


	// *********************************************
	// Quantize the noisy image
	// *********************************************
	for (int j = 0; j < lena.rows; j++) {
		for (int i = 0; i < lena.cols; i++) {

			int grayVal = (int) noisy.at<uchar>(i, j);

			int mod32 = grayVal % 32;
			int difference = 16 - mod32;
			int newVal = grayVal + difference;

			// assign quantized value to quantized image
			noisyQuan.at<uchar>(i, j) = (uchar) newVal;

			// print
			// cout << "noisy(" << i << ", " << j << "): " << (int) noisy.at<uchar>(i, j);
			// cout << "\tnoisyQuan(" << i << "," << j << "): " << (int) noisyQuan.at<uchar>(i, j) << endl;
			// cout << endl;
		}
	}

	cout << "*** Mean Squared Error between noisy image and quantized noisy image ***" << endl;
	findMSE(noisy, noisyQuan);

	// show
	imshow("Noisy", noisy);
	imshow("Noisy Quantized", noisyQuan);
	while(waitKey(30) != 27);
	destroyAllWindows();

	// write
	imwrite("../images(#7)/noisyQuantized.png", noisyQuan);


	// *********************************************
	// Apply the lowpass filter kernel
	// *********************************************

	float values[9] = {1, 2, 1, 2, 4, 2, 1, 2, 1};

	Mat kernel(3, 3, CV_32F, values);
	float scal = 1.0f / 16.0f;
	kernel *= scal;

	// filter - (-1 for default depth)
	filter2D(quan, quanFiltered, -1, kernel);

	// filter noisy and noisy quantized images
	filter2D(noisy, noisyFiltered, -1, kernel);
	filter2D(noisyQuan, noisyQuanFiltered, -1, kernel);

	cout << "*** Mean Squared Error between original image and quantized/filtered image ***" << endl;
	findMSE(lena, quanFiltered);

	cout << "*** Mean Squared Error between original image and noisy/filtered image ***" << endl;
	findMSE(lena, noisyFiltered);

	cout << "*** Mean Squared Error between original image and noisy/quantized/filtered image ***" << endl;
	findMSE(lena, noisyQuanFiltered);



	imshow("Quantized", quan);
	imshow("Quantized - Filtered", quanFiltered);
	imshow("Noisy - Filtered", noisyFiltered);
	imshow("Noisy/Quantized - Filtered", noisyQuanFiltered);
	while(waitKey(30) != 27);
	destroyAllWindows();

	// write
	imwrite("../images(#7)/quantizedFiltered.png", quanFiltered);
	imwrite("../images(#7)/noisyFiltered.png", noisyFiltered);
	imwrite("../images(#7)/noisyQuantizedFiltered.png", noisyQuanFiltered);


	return 0;
}



int main() {

	cout << endl;

	// problemSix();

	//problemSeven();

	laplacian();

	// highPassEqualization();

	return 0;
}

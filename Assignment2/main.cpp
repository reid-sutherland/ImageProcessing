#include <iostream>
#include <ctime>
#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <opencv2/xfeatures2d/nonfree.hpp>

using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

void QuestionOne() {
    // timer stuff
    clock_t startTime, checkTime, timePassed;
    double millis_sift, millis_surf;


    Mat img1 = imread("Assig2_Supp/Q1-Invariant-Testing/box.png", IMREAD_GRAYSCALE);
    Mat img2 = imread("Assig2_Supp/Q1-Invariant-Testing/box_in_scene.png", IMREAD_GRAYSCALE);

    if (img1.empty()) {
        cerr << "Error: Image 1 empty" << endl;
        return;
    }
    if (img2.empty()) {
        cerr << "Error: Image 2 empty" << endl;
        return;
    }

    vector<KeyPoint> keypoints1, keypoints2;
    Mat descriptors1, descriptors2;
    vector<DMatch> matches1, matches2;
    Ptr<DescriptorMatcher> matcher = BFMatcher::create("BruteForce");


    // SIFT detection
    // start timer
    startTime = clock();

    // create pointers
    Ptr<SIFT> sift1 = SIFT::create();
    Ptr<SIFT> sift2 = SIFT::create();

    // start detection
    sift1->detect(img1, keypoints1);
    sift2->detect(img2, keypoints2);

    // compute descriptors
    sift1->compute(img1, keypoints1, descriptors1);
    sift2->compute(img2, keypoints2, descriptors2);

    // find matches
    matcher->match(descriptors1, matches1);
    matcher->match(descriptors2, matches2);

    // draw keypoints
    Mat img1_sift, img2_sift;
    drawKeypoints(img1, keypoints1, img1_sift);
    drawKeypoints(img2, keypoints2, img2_sift);

    // draw matches
    Mat sift_matches;
    drawMatches(img1, keypoints1, img2, keypoints2)
    

    // stop timer, record time
    checkTime = clock();
    timePassed = checkTime - startTime;
    millis_sift = timePassed / (double)CLOCKS_PER_SEC * 1000;


    // SURF detection
    // start timer
    startTime = clock();

    // create pointers
    int hessianThresh = 300;
    Ptr<SURF> surf1 = SURF::create(hessianThresh);
    Ptr<SURF> surf2 = SURF::create(hessianThresh);

    // start detection
    surf1->detect(img1, keypoints1);
    surf2->detect(img2, keypoints2);

    // draw keypoints
    Mat img1_surf, img2_surf;
    drawKeypoints(img1, keypoints1, img1_surf);
    drawKeypoints(img2, keypoints2, img2_surf);

    // stop timer, record time
    checkTime = clock();
    timePassed = checkTime - startTime;
    millis_surf = timePassed / (double)CLOCKS_PER_SEC * 1000;


    // show results
    cout << "SIFT took approx. " << millis_sift << " milliseconds." << endl;
    cout << "SURF took approx. " << millis_surf << " milliseconds." << endl;
    imshow("Image 1 SIFT", img1_sift);
    imshow("Image 2 SIFT", img2_sift);
    imshow("Image 1 SURF", img1_surf);
    imshow("Image 2 SURF", img2_surf);

    waitKey(0);
}

int main() {

    QuestionOne();

    return 0;
}
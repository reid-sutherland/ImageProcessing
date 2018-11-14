#include <iostream>
#include <ctime>
#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <opencv2/xfeatures2d/nonfree.hpp>

using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

// limit matches by an arbitrary distance
// store them all into one vector
vector<DMatch> matchCeiling(vector<vector<DMatch>> in_matches, float ceiling) {
    vector<DMatch> goodMatches;

    clock_t startTime = clock();

    // loop through all matches
    for (int i = 0; i < in_matches.size(); i++) {
        for (int j = 0; j < in_matches[i].size(); j++) {
            if (in_matches[i][j].distance <= ceiling)
                goodMatches.push_back(in_matches[i][j]);
        }
    }

    clock_t timePassed = (clock() - startTime);
    double millis = timePassed / (double)CLOCKS_PER_SEC * 1000;
    cout << "matchCeiling time:  " << millis << "ms" << endl;
}

// limit matches by an arbitrary distance
// preserve grouping by each query descriptor
vector<vector<DMatch>> matchCeiling2(vector<vector<DMatch>> in_matches, float ceiling) {
    vector<vector<DMatch>> goodMatches;

    clock_t startTime = clock();
    // loop through all matches
    for (int i = 0; i < in_matches.size(); i++) {
        vector<DMatch> row;
        for (int j = 0; j < in_matches[i].size(); j++) {
            if (in_matches[i][j].distance <= ceiling)
                row.push_back(in_matches[i][j]);
        }
        goodMatches.push_back(row);
    }

    clock_t timePassed = (clock() - startTime);
    double millis = timePassed / (double)CLOCKS_PER_SEC * 1000;
    cout << "matchCeiling2 time:  " << millis << "ms" << endl;
}

// find best(ish) matches from a 2d array of matches
vector<DMatch> findBestMatches(vector<vector<DMatch>> in_matches, int percentile) {
    // container for good matches
    vector<DMatch> goodMatches;
    // insert the min value (first element) of first row of DMatches
    goodMatches.push_back(in_matches[0][0]);
    float thresh = goodMatches[0].distance;     // threshold for good matches
    int target = (int) in_matches.size() / (100/percentile);    // target number of good matches
    // mask of zeros, set to non-zero when match found
    Mat mask = Mat::zeros((int)in_matches.size(), (int)in_matches[0].size(), CV_32S);

    // get best n matches, put into list
    while (goodMatches.size() < target) {
//        cout << "\nLooking for good matches..." << endl;
//        cout << "Threshold: " << thresh << endl;
        float newThreshDist = 0.0;
        for (int i = 0; i < in_matches.size(); i++) {
            for (int j = 0; j < in_matches[i].size(); j++) {
                float dist = in_matches[i][j].distance;
                // if near (but less than) threshold, set new anchor for next run (if needed)
                if (dist - thresh < 10 && dist - thresh > 0 && newThreshDist == 0.0) { //don't if new thresh has been found
                    newThreshDist = dist;
                }
                // if already good match, don't add
                if (mask.at<int>(i, j) != 0)
                    continue;
                // otherwise, add the good match, mark it on the mask
                if (dist <= thresh) {
                    goodMatches.push_back(in_matches[i][j]);
                    mask.at<int>(i, j) = 1;
                }

            }
        }
//        cout << "\nNum good matches found: " << goodMatches.size() << endl;
//        cout << "Target: " << target << endl;
        thresh = newThreshDist;
    }

    // sort list of best matches and limit to preferred size (insertion sort)
    int i = 1;
    while (i < goodMatches.size()) {    // alt insertion sort, big values first
        DMatch temp = goodMatches[i];
        int j = i-1;
        while (j >= 0 && goodMatches[j].distance > temp.distance) {
            goodMatches[j+1] = goodMatches[j];
            j--;
        }
        goodMatches[j+1] = temp;
        i++;
    }

    // resize to fit target
    goodMatches.resize(target);

    return goodMatches;
}

// use insertion sort to sort a list of matches
vector<DMatch> sortMatches(vector<DMatch> in_matches, int percentile) {
    int i = 1;
    while (i < in_matches.size()) {
        DMatch temp = in_matches[i];
        int j = i-1;
        while (j >= 0 && in_matches[j].distance > temp.distance) {
            in_matches[j+1] = in_matches[j];
            j--;
        }
        in_matches[j+1] = temp;
        i++;
    }

    // resize by percentile
    unsigned long target = in_matches.size() / (100/percentile);
    in_matches.resize(target);

    return in_matches;
}

void QuestionOne() {
    //! Load images
    Mat img1, img2, img3;
    int percentile;     // percentile of matches to use as good matches
    int hessianThresh;    // hessian threshold for surf detection

    // 0  - box
    // 1  - DSC_1
    // 12 - DSC_2
    // 2  - iphone_1
    // 22 - iphone_2
    // 3  - label
    int subject = 3;
    string subjectName;

    switch (subject) {
        // box
        case 0:
            img1 = imread("Assig2_Supp/Q1-Invariant-Testing/box.png", IMREAD_GRAYSCALE);
            img2 = imread("Assig2_Supp/Q1-Invariant-Testing/box_in_scene.png", IMREAD_GRAYSCALE);
            subjectName = "box";
            percentile = 5;
            hessianThresh = 500;
            break;

        // DSC
        case 1:
            img1 = imread("Assig2_Supp/Q1-Invariant-Testing/DSC_0047.JPG", IMREAD_GRAYSCALE);
            img2 = imread("Assig2_Supp/Q1-Invariant-Testing/DSC_0053.JPG", IMREAD_GRAYSCALE);
            img2 = imread("Assig2_Supp/Q1-Invariant-Testing/DSC_0057.JPG", IMREAD_GRAYSCALE);
            subjectName = "DSC1";
            percentile = 5;
            hessianThresh = 500;
            break;

        // DSC_2
        case 12:
            img1 = imread("Assig2_Supp/Q1-Invariant-Testing/DSC_0047.JPG", IMREAD_GRAYSCALE);
            img2 = imread("Assig2_Supp/Q1-Invariant-Testing/DSC_0057.JPG", IMREAD_GRAYSCALE);
            subjectName = "DSC2";
            percentile = 5;
            hessianThresh = 500;
            break;

        // iphone
        case 2:
            img1 = imread("Assig2_Supp/Q1-Invariant-Testing/iphone_object.jpg", IMREAD_GRAYSCALE);
            img2 = imread("Assig2_Supp/Q1-Invariant-Testing/iphone1.jpg", IMREAD_GRAYSCALE);
            subjectName = "iphone1";
            percentile = 7;
            hessianThresh = 500;
            break;

        // iphone_2
        case 22:
            img1 = imread("Assig2_Supp/Q1-Invariant-Testing/iphone_object.jpg", IMREAD_GRAYSCALE);
            img2 = imread("Assig2_Supp/Q1-Invariant-Testing/iphone2.jpg", IMREAD_GRAYSCALE);
            subjectName = "iphone2";
            percentile = 7;
            hessianThresh = 500;
            break;

        // label
        case 3:
            img1 = imread("Assig2_Supp/Q1-Invariant-Testing/label_object.jpg", IMREAD_GRAYSCALE);
            img2 = imread("Assig2_Supp/Q1-Invariant-Testing/label.jpg", IMREAD_GRAYSCALE);
            subjectName = "label";
            percentile = 15;
            hessianThresh = 500;
            break;

        default:
            cout << "Error: Need to specify a subject." << endl;
            return;
    }


    if (img1.empty()) {
        cerr << "Error: Image 1 empty" << endl;
        return;
    }
    if (img2.empty()) {
        cerr << "Error: Image 2 empty" << endl;
        return;
    }

    //! Setup
    // Timer stuff
    clock_t startTime, checkTime, timePassed;
    double millis_sift1, millis_surf1;
//    double millis_sift2, millis_surf2;

    // Containers
    vector<KeyPoint> keypoints1, keypoints2;
    Mat descriptors1, descriptors2;
    vector<DMatch> matches;
    vector<vector<DMatch>> matches2;
	BFMatcher matcher(NORM_L2);
	float totalDistance, avgDistance;

	// Images
	Mat img1_sift, img2_sift, img_matches_sift;
	Mat img1_surf, img2_surf, img_matches_surf;
//    Mat img1_sift2, img2_sift2, img_matches_sift2;
//    Mat img1_surf2, img2_surf2, img_matches_surf2;


    //! SIFT detection
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

    // find matches across combined two images
    matcher.match(descriptors1, descriptors2, matches);
//    matcher.radiusMatch(descriptors1, descriptors2, matches2, img1.rows+img1.cols+img2.rows+img2.cols);

    // reduce list of matches to best 5% only (% based on num keypoints in object image
    // sort list of matches, best first (insertion sort)
    vector<DMatch> goodMatches;
    goodMatches = sortMatches(matches, percentile);
//    goodMatches = findBestMatches(matches2, percentile);

    // draw keypoints
    drawKeypoints(img1, keypoints1, img1_sift);
    drawKeypoints(img2, keypoints2, img2_sift);

    // draw matches
    drawMatches(img1, keypoints1, img2, keypoints2, goodMatches, img_matches_sift);

    // stop timer, record time
    checkTime = clock();
    timePassed = checkTime - startTime;
    millis_sift1 = timePassed / (double)CLOCKS_PER_SEC * 1000;

    // calculate distance
    totalDistance = 0;
    avgDistance = 0;
    for (int i = 0; i < goodMatches.size(); i++) {
        totalDistance += goodMatches[i].distance;
    }
    avgDistance = totalDistance / goodMatches.size();

    // print
    cout << endl;
    cout << "******** " << subjectName << " ********" << endl << endl;
    cout << "_____SIFT_____" << endl;
    cout << "Quantity" << endl;
    cout << " Object features:  " << keypoints1.size() << endl;
    cout << " _in_space features:  " << keypoints2.size() << endl;
    cout << " Matches:  " << matches.size() << endl;
    cout << " Good matches:  " << goodMatches.size() << endl;
    cout << "Performance" << endl;
    cout << " Computation time:  " << millis_sift1 << "ms" << endl;
    cout << " Matching accuracy (total distance):  " << totalDistance << endl;
    cout << " Matching accuracy (average distance):  " << avgDistance << " per match" << endl;
    cout << endl;



    //! SURF detection
    // start timer
    startTime = clock();

    // create pointers
    Ptr<SURF> surf1 = SURF::create(hessianThresh);
    Ptr<SURF> surf2 = SURF::create(hessianThresh);

    // start detection
    surf1->detect(img1, keypoints1);
    surf2->detect(img2, keypoints2);

    // compute matches
    surf1->compute(img1, keypoints1, descriptors1);
    surf2->compute(img2, keypoints2, descriptors2);

    // find matches
    matcher.match(descriptors1, descriptors2, matches);
    // reduce list of matches to best x% only (% based on num keypoints in object image)
    // sort list of matches, best first (insertion sort)
    goodMatches = sortMatches(matches, percentile);
//    goodMatches = findBestMatches(matches2, percentile);

    // draw keypoints
    drawKeypoints(img1, keypoints1, img1_surf);
    drawKeypoints(img2, keypoints2, img2_surf);

    // draw matches
    drawMatches(img1, keypoints1, img2, keypoints2, goodMatches, img_matches_surf);

    // stop timer, record time
    checkTime = clock();
    timePassed = checkTime - startTime;
    millis_surf1 = timePassed / (double)CLOCKS_PER_SEC * 1000;

    // calculate distance
    totalDistance = 0;
    avgDistance = 0;
    for (int i = 0; i < goodMatches.size(); i++) {
        // SURF distance is normalized, so multiply by total number of matches
        totalDistance += (goodMatches[i].distance * matches.size());
    }
    avgDistance = totalDistance / goodMatches.size();

    // SURF distance is normalized, so multiply by number of matches


    // print
    cout << "_____SURF_____" << endl;
    cout << "Quantity" << endl;
    cout << " Object features:  " << keypoints1.size() << endl;
    cout << " _in_space features:  " << keypoints2.size() << endl;
    cout << " Matches:  " << matches.size() << endl;
    cout << " Good matches:  " << goodMatches.size() << endl;
    cout << "Performance" << endl;
    cout << " Computation time:  " << millis_surf1 << "ms" << endl;
    cout << " Matching accuracy (total distance):  " << totalDistance << endl;
    cout << " Matching accuracy (average distance):  " << avgDistance << " per match" << endl;


    //! Do it all again if there is a second _in_scene image
//    if (!img3.empty()) {
//        //! SIFT detection
//        // start timer
//        startTime = clock();
//        // create pointers
//        Ptr<SIFT> sift3 = SIFT::create();
//        Ptr<SIFT> sift4 = SIFT::create();
//        // start detection
//        sift3->detect(img1, keypoints1);
//        sift4->detect(img3, keypoints2);
//        // compute descriptors
//        sift3->compute(img1, keypoints1, descriptors1);
//        sift4->compute(img3, keypoints2, descriptors2);
//        // find matches across combined two images
//        matcher.match(descriptors1, descriptors2, matches);
//        matcher.radiusMatch(descriptors1, descriptors2, matches2, img1.rows+img1.cols+img3.rows+img3.cols);
//        // reduce list of matches to best 5% only (% based on num keypoints in object image)
//        // use same percentile as before
//        // sort list of matches, best first (insertion sort)
//        goodMatches = sortMatches(matches, percentile);
////    goodMatches = findBestMatches(matches2, percentile);
//        // draw keypoints
//        drawKeypoints(img1, keypoints1, img1_sift2);
//        drawKeypoints(img3, keypoints2, img2_sift2);
//        // draw matches
//        drawMatches(img1, keypoints1, img3, keypoints2, goodMatches, img_matches_sift2);
//        // stop timer, record time
//        checkTime = clock();
//        timePassed = checkTime - startTime;
//        millis_sift2 = timePassed / (double)CLOCKS_PER_SEC * 1000;
//
//        //! SURF detection
//        // start timer
//        startTime = clock();
//        // create pointers
//        Ptr<SURF> surf3 = SURF::create(hessianThresh);  // use same Hessian Threshold as before
//        Ptr<SURF> surf4 = SURF::create(hessianThresh);
//        // start detection
//        surf3->detect(img1, keypoints1);
//        surf4->detect(img3, keypoints2);
//        // compute matches
//        surf3->compute(img1, keypoints1, descriptors1);
//        surf4->compute(img3, keypoints2, descriptors2);
//        // find matches
//        matcher.match(descriptors1, descriptors2, matches);
//        // reduce list of matches to best 5% only (% based on num keypoints in object image)
//        // sort list of matches, best first (insertion sort)
//        goodMatches = sortMatches(matches, percentile);
////    goodMatches = findBestMatches(matches2, percentile);
//        // draw keypoints
//        drawKeypoints(img1, keypoints1, img1_surf2);
//        drawKeypoints(img3, keypoints2, img2_surf2);
//        // draw matches
//        drawMatches(img1, keypoints1, img3, keypoints2, goodMatches, img_matches_surf2);
//        // stop timer, record time
//        checkTime = clock();
//        timePassed = checkTime - startTime;
//        millis_surf2 = timePassed / (double)CLOCKS_PER_SEC * 1000;
//    }


    // show results
    imshow("Object Keypoints - SIFT", img1_sift);
    imshow("_in_scene Keypoints - SIFT", img2_sift);
    imshow("Matches SIFT", img_matches_sift);
    imshow("Object Keypoints - SURF", img1_surf);
    imshow("_in_scene Keypoints - SURF", img2_surf);
    imshow("Matches SURF", img_matches_surf);

    waitKey(0);

    // create output names
    string sift_name1, sift_name2, sift_matches, surf_name1, surf_name2, surf_matches;
    sift_name1 = "output/";
    sift_name1 += subjectName;
    sift_name1 += "/";
    sift_name2 = sift_name1;  sift_matches = sift_name1;  surf_name1 = sift_name1;  surf_name2 = sift_name1;  surf_matches = sift_name1;

    sift_name1 += "SIFT";
    sift_name2 += "SIFT";
    sift_matches += "SIFT";
    surf_name1 += "SURF";
    surf_name2 += "SURF";
    surf_matches += "SURF";

    sift_name1 += "_object.png";
    surf_name1 += "_object.png";
    sift_name2 += "_in_space.png";
    surf_name2 += "_in_space.png";
    sift_matches += "_matches.png";
    surf_matches += "_matches.png";

    // write images
    imwrite(sift_name1, img1_sift);
    imwrite(sift_name2, img2_sift);
    imwrite(sift_matches, img_matches_sift);
    imwrite(surf_name1, img1_surf);
    imwrite(surf_name2, img2_surf);
    imwrite(surf_matches, img_matches_surf);

    cout << sift_name1 << endl;
    cout << surf_name1 << endl;


//    if (!img3.empty()) {
//        cout << "\nSecond _in_scene image:  " << endl;
//        cout << "SIFT took approx. " << millis_sift2 << " milliseconds." << endl;
//        cout << "SURF took approx. " << millis_surf2 << " milliseconds." << endl;
//        imshow("Object Keypoints - SIFT", img1_sift2);
//        imshow("_in_scene Keypoints - SIFT", img2_sift2);
//        imshow("Matches SIFT", img_matches_sift2);
//        imshow("Object Keypoints - SURF", img1_surf2);
//        imshow("_in_scene Keypoints - SURF", img2_surf2);
//        imshow("Matches SURF", img_matches_surf2);
//
//        waitKey(0);
//    }
}

int main() {

//    clock_t sortTime = clock();
//    cout << "Sort on matches took approx. " << (clock() - sortTime) / (double)CLOCKS_PER_SEC*1000 << " milliseconds." << endl;

    QuestionOne();

    return 0;
}
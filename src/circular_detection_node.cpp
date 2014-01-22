/*
 * circular_detection_node.cpp
 *
 *  Created on: Jan 22, 2014
 *      Author: windsdon
 */

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>

#include <vector>

using namespace cv;
using namespace std;

const float PI = 3.141592653;

#define RAYS 8
#define MAX_ERROR 0.1
#define MAX_HUE 180

Mat *imageDraw = NULL;

bool within(Vec3b &val, int *lower, int *upper, int channels, int shift) {
	for (int i = 0; i < channels; i++) {
		uchar v;
		if (!i) {
			v = (val[i] + shift) % MAX_HUE;
		} else {
			v = val[i];
		}
		if (v < lower[i] || v > upper[i]) {
			return false;
		}
	}
	return true;
}

typedef Point2f Vector2f;

class MathHelper {
	public:
		static float size2(Vector2f& v) {
			return (pow(v.x, 2) + pow(v.y, 2));
		}

		static float size(Vector2f &v) {
			return sqrt(size2(v));
		}
		static Vector2f perp(Vector2f& v) {
			return Vector2f(-v.y, v.x);
		}

		static Point2f circleCentre(vector<Point2f> &points) {
			int pcount = points.size();

			float cx = 0;
			float cy = 0;
			int interCount = 0;

			for (int i = 0; i < pcount; i++) {
				Point2f b = points[(i + 1) % pcount];
				Point2f a = points[i];
				Vector2f r = b - a;
				Point2f p0 = a + r * 0.5;
				Point2f p1 = p0 + perp(r);

				for (int j = 0; j < pcount; j++) {
					if (i == j) {
						continue;
					}

					Point2f b2 = points[(j + 1) % pcount];
					Point2f a2 = points[j];
					Vector2f r2 = b2 - a2;
					Point2f q0 = a2 + r2 * 0.5;
					Point2f q1 = q0 + perp(r2);

					Vector2f intersect;
					if (lineIntersect(p0, p1, q0, q1, intersect)) {
						cx += intersect.x;
						cy += intersect.y;

						interCount++;
					}

				}
			}

			if (interCount > 0) {
				return Point2f(cx / interCount, cy / interCount);
			} else {
				return Point2f();
			}

		}

		static float errorPoint(Point2f &centre, vector<Point2f> &points) {
			float e = 0;
			for (vector<Point2f>::iterator it = points.begin(); it != points.end(); ++it) {
				Vector2f v = (*it) - centre;
				e += size(v);
			}

			return (e / points.size());
		}

		static float radius(Point2f &centre, vector<Point2f> &points, float &error) {
			float r = 0;

			vector<float> rc;

			for (vector<Point2f>::iterator it = points.begin(); it != points.end(); ++it) {
				Vector2f ri = *it - centre;

				float rl = size(ri);

				rc.push_back(rl);

				r += rl;
			}

			r /= points.size();

			float dev = 0;

			for (vector<float>::iterator it = rc.begin(); it != rc.end(); ++it) {
				dev += abs(*it - r);
			}

			error = dev / rc.size();

			return r;
		}

		static bool lineIntersect(Point2f &p0, Point2f &p1, Point2f &q0, Point2f &q1, Point2f& intersect) {
			Vector2f u = p1 - p0;
			Vector2f v = q1 - q0;

			float det = (u.x * v.y - u.y * v.x);

			if (det == 0) {
				return false;
			}

			float s = ((p0.y - q0.y) * v.x + (q0.x - p0.x) * v.y) / det;
			float t = ((q0.x - p0.x) * u.y + (p0.y - q0.y) * u.x) / det;

			if (s < 0 || t < 0) {
				//special case for the circle centre
				return false;
			}

			intersect = p0 + s * u;

			return true;
		}
};

Point2f getEdge(Mat &image, int *lower, int *upper, int sx, int sy, float angle, int shift) {
	int x, y;
	int px = sx, py = sy;
	for (int r = 0;; r++) {
		x = sx + r * cos(angle);
		y = sy + r * sin(angle);

		if (x < 0 || y < 0 || x >= image.cols || y >= image.rows) {
			break;
		}

		Vec3b val = image.at<Vec3b>(y, x);
		if (!within(val, lower, upper, 3, shift)) {
			break;
		}

		px = x;
		py = y;
	}

	return Point2f(px, py);
}
vector<Point2f> cache;
float getCircle(Mat &image, int *lower, int *upper, float sx, float sy, float *cx, float *cy, float *r, int shift, int rays) {
	float angleStep = 2 * PI / rays;

	vector<Point2f> edges;

	cache.clear();

	for (int step = 0; step < rays; step++) {
		float angle = angleStep * step;

		edges.push_back(getEdge(image, lower, upper, sx, sy, angle, shift));

		cache.push_back(edges.back());

		//circle(*imageDraw, edges.back(), 1, Scalar(0, 255, 255), -1);
	}

	Point2f centre = MathHelper::circleCentre(edges);

	float error;

	*r = MathHelper::radius(centre, edges, error);
	*cx = centre.x;
	*cy = centre.y;

	return error;
}

bool findCircle(Mat &image, int *lower, int *upper, float *cx, float *cy, float *r, float minR, int boxStartX, int boxStartY, int boxEndX, int boxEndY, int shift) {
	//this function trusts that the box is contained in the image, and performs no checks

	int width = boxEndX - boxStartX;
	int height = boxEndY - boxStartY;

	float p = 1 / 4.0;

	int xstep, ystep;

	do {
		xstep = ceil(width * p);
		ystep = ceil(height * p);

		for (int x = boxStartX; x <= boxEndX; x += xstep) {
			for (int y = boxStartY; y <= boxEndY; y += ystep) {
				Vec3b pixel = image.at<Vec3b>(y, x);

				if (within(pixel, lower, upper, 3, shift)) {
					float cx2, cy2, r2;

					float error = getCircle(image, lower, upper, x, y, &cx2, &cy2, &r2, shift, RAYS);

					cout << "r: " << r2 << " \033[31mError: " << error << ", " << (error / r2) << "\033[0m" << endl;

					if (error / r2 < MAX_ERROR && r2 >= minR) {
						*r = r2;
						*cx = cx2;
						*cy = cy2;

						for(vector<Point2f>::iterator it = cache.begin(); it != cache.end(); ++it){
							circle(*imageDraw, *it, 2, Scalar(0, 255, 255), -1);
						}

						//circle(*imageDraw, Point2i(x, y), 2, Scalar(0, 255, 0), -1);

						return true;
					}
				}
			}
		}

		p /= 2;
	} while (xstep > 1 || ystep > 1);

	return false;
}

Mat *imageRaw = NULL;

void onMouse(int event, int x, int y, int, void*) {
	if (event != EVENT_LBUTTONDOWN) {
		return;
	}

	Vec3b pixel = imageRaw->at<Vec3b>(y, x);

	cout << pixel << endl;

}

void makeBounds(int *target, int *delta, int *upper, int *lower, int &hueShift) {
	hueShift = 0;

	int low, high;
	for (int i = 0; i < 3; i++) {
		high = upper[i] = target[i] + delta[i];
		low = lower[i] = target[i] - delta[i];

		if (!i) {
			if (low < 0 && abs(hueShift) < -low) {
				hueShift = -low;
			}

			if (high >= MAX_HUE && (high - MAX_HUE) > abs(hueShift)) {
				hueShift = MAX_HUE - high;
			}

			upper[i] += hueShift;
			lower[i] += hueShift;

			if (upper[i] >= MAX_HUE) {
				upper[i] = MAX_HUE - 1;
			}
		}

		if(upper[i] > 255){
			upper[i] = 255;
		}

		if (lower[i] < 0) {
			lower[i] = 0;
		}

	}
}

int main(int argc, char **argv) {

	Mat image = imread("image.png");

	medianBlur(image, image, 13);

	imageDraw = new Mat();
	imageRaw = &image;

	image.copyTo(*imageDraw);

	cvtColor(image, image, CV_BGR2HSV);

	int upper[3];
	int lower[3];
	int delta[] = {10, 20, 20};
	int target[] = {13, 214, 251};
	int shift;

	makeBounds(target, delta, upper, lower, shift);

	Vec3b low(lower[0], lower[1], lower[2]);
	Vec3b up(upper[0], upper[1], upper[2]);

	cout << low << ", " << up << ", " << shift << endl;

	float cx, cy, r;

	if (findCircle(image, lower, upper, &cx, &cy, &r, 40, 0, 0, image.cols - 1, image.rows - 1, shift)) {
		cout << "Found circle" << endl;
		circle(*imageDraw, Point2f(cx, cy), 4, Scalar(0, 0, 255), -1);

		circle(*imageDraw, Point2f(cx, cy), r, Scalar(0, 0, 255), 2);
	} else {
		cout << "Didn't find circle" << endl;
	}

	imshow("Image", *imageDraw);

	//waitKey(0);

	setMouseCallback("Image", &onMouse);

	while (true) {
		if (waitKey(3) == 'x') {
			break;
		}
	}

}


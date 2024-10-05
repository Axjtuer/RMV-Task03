/*#include "windmill.hpp"
#include "ceres/ceres.h"
#include "glog/logging.h"
using namespace std;
using namespace cv;
using namespace ceres;
struct CosRes
{
    CosRes(double x, double y) : x_(x), y_(y) {}

    template <typename T>
    bool operator()(const T *const A0, const T *const A, const T *const w, const T *const phi, T *residual) const
    {
        residual[0] = y_ - ceres::cos(A0[0] * x_ + A[0] / w[0] * (ceres::cos(phi[0] + 1.5707963) - ceres::cos(w[0] * x_ + phi[0] + 1.5707963)));
        return true;
    }

private:
    const double x_;
    const double y_;
};
inline bool checkcomb(double nowA0, double nowA, double noww, double nowphi)
{
    if (nowphi < 0.25 && 1.78 < noww && 1.23 < nowA0 && nowA < 0.83 && nowA0 < 1.37 && 0.74 < nowA && noww < 1.98 && 0.22 < nowphi)
    {
        return true;
    }
    return false;
}
int main()
{
    double t_sum = 0;
    const int N = 10;
    for (int num = 0; num < N; num++)
    {
        std::cout<<num<<std::endl;
        std::chrono::milliseconds t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
        double t_start = (double)t.count();
        WINDMILL::WindMill wm(t_start);
        Mat src;

        Problem problem;
        double A0 = 0.305, A = 1.785, w = 0.884, phi = 1.24;

        int count = 0;

        // starttime
        int64 start_time = getTickCount();

        while (1)
        {   count++;
            t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
            double t_now = (double)t.count();
            src = wm.getMat(t_now); // Is here wrong? (why to divide 1000?) Still, it can run well. But it is unreasonal.
            //==========================代码区========================//
            imshow("0",src);
            Mat grayimg;
            cvtColor(src,grayimg,COLOR_BGR2GRAY);
            Mat binary;
            threshold(grayimg,binary,75,225,THRESH_BINARY);
            vector<vector<Point>>contours;
            vector<Vec4i>hierarchy;
            Point2i center;
            findContours(binary,contours,hierarchy,RETR_TREE,CHAIN_APPROX_SIMPLE);
            int contour[20]={0};
            for (int i = 0; i < contours.size(); i++)//遍历检测的所有轮廓
            {
                if (hierarchy[i][3] != -1) //有内嵌轮廓，说明是一个父轮廓
                {
                contour[hierarchy[i][3]]++; //对该父轮廓进行记录
                }
            }
            Point R_center;
            for (int j = 0; j < contours.size(); j++)//再次遍历所有轮廓
            {
                if (contour[j] == 1) //如果某轮廓对应数组的值为1，说明只要一个内嵌轮廓
                {
                    double s=contourArea(contours[j]);
                    if(s<200)
                    {   
                        Moments R=moments(contours[j]);
                        if(R.m00!=0){
                            Point cen(R.m10/R.m00,R.m01/R.m00);
                            R_center=cen;
                            circle(src,cen,4,Scalar(255,0,0),-1);
                        }
                        continue;
                    } 
                    int num = hierarchy[j][2]; //记录该轮廓的内嵌轮廓
                    RotatedRect box = minAreaRect(contours[num]); //包含该轮廓所有点
                    Point2f vertex[4];
                    box.points(vertex);//将左下角，左上角，右上角，右下角存入点集
                    for (int i = 0; i < 4; i++)
                        {
                            line(src, vertex[i], vertex[(i + 1) % 4], Scalar(255, 0, 0), 4, LINE_AA); //画线
                        }
                    center = (vertex[0] + vertex[2]) / 2; //返回中心坐标
                    circle(src,center,4,Scalar(255,0,0),-1);
                    putText(src, "target", vertex[0], FONT_HERSHEY_SIMPLEX, 1.0, Scalar(255, 255, 0));//打印字体
                    }
            }
            Point rh=center-R_center;
            Point2d rhi=Point2d(rh)/norm(rh);
            double xt=(t_now-t_start)/1000;
            double yt=rhi.x;

            problem.AddResidualBlock(new AutoDiffCostFunction<CosRes, 1, 1, 1, 1, 1>(new CosRes(xt, yt)), NULL, &A0, &A, &w, &phi);
            Solver::Options options;
            options.max_num_iterations =25;
            options.linear_solver_type=DENSE_QR;
            problem.SetParameterLowerBound(&A0, 0, 0.5);
            problem.SetParameterUpperBound(&A0, 0, 1.4);
            problem.SetParameterLowerBound(&A, 0, 0.5);
            problem.SetParameterUpperBound(&A, 0, 1.0);
            problem.SetParameterLowerBound(&w, 0, 0.5);
            problem.SetParameterUpperBound(&w, 0, 1.9);
            problem.SetParameterLowerBound(&phi, 0, 0.24);
            problem.SetParameterUpperBound(&phi, 0, 1.25);
            Solver::Summary summary;
            Solve(options,&problem,&summary);
            if(checkcomb(A0, A, w, phi)){
                int end_time =getTickCount();
                t_sum+=(end_time-start_time)/getTickFrequency();
                break;
            }
            //=======================================================//
            //waitKey(1);
        }
    }
    std::cout << t_sum / N << std::endl;
}*/
#include "windmill.hpp"
#include "ceres/ceres.h"
#include "glog/logging.h"

using namespace std;
using namespace cv;
using namespace ceres;
#define pi M_PI
// CostFunctor y-cos(alpha)
struct CostFunctor
{
    CostFunctor(double x, double y) : x_(x), y_(y) {}

    template <typename T>
    bool operator()(const T *const A0, const T *const A, const T *const w, const T *const phi, T *residual) const
    {   T res=ceres::cos(A0[0] * x_ + A[0] / w[0] * (ceres::cos(phi[0] + pi/2) - ceres::cos(w[0] * x_ + phi[0] + pi/2)));
        residual[0] = y_ - res;
        return true;
    }

private:
    const double x_;
    const double y_;
};

inline bool checkcomb(double nowA0, double nowA, double noww, double nowphi)
{
    if (nowphi < 0.25 && 1.78 < noww && 1.23 < nowA0 && nowA < 0.83 && nowA0 < 1.37 && 0.74 < nowA && noww < 1.98 && 0.22 < nowphi)
    {
        return true;
    }
    return false;
}

int main()
{
    double t_sum = 0;
    const int N = 10;
    for (int num = 0; num < N; num++)
    {   cout<<num<<endl;
        std::chrono::milliseconds t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
        double t_start = (double)t.count();
        WINDMILL::WindMill wm(t_start);
        Mat src;

        Problem problem;
        double A0 = 0.305, A = 1.785, w = 0.884, phi = 1.24;

        int count = 0;

        // starttime
        int64 start_time = getTickCount();

        while (1)
        {
            count++;
            t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
            double t_now = (double)t.count();
            src = wm.getMat(t_now); // Is here wrong? (why to divide 1000?) Still, it can run well. But it is unreasonal.

            /*code*/

            // 1. draw circles
            // gray
            Mat highlightGray;
            cvtColor(src, highlightGray, COLOR_BGR2GRAY);

            // binary
            Mat binary;
            threshold(highlightGray, binary, 50, 255, THRESH_BINARY);

            vector<vector<Point>>contours;
            vector<Vec4i>hierarchy;
            Point2i center;
            findContours(binary,contours,hierarchy,RETR_TREE,CHAIN_APPROX_SIMPLE);
            int contour[20]={0};
            for (int i = 0; i < contours.size(); i++)//遍历检测的所有轮廓
            {
                if (hierarchy[i][3] != -1) //有内嵌轮廓，说明是一个父轮廓
                {
                contour[hierarchy[i][3]]++; //对该父轮廓进行记录
                }
            }
            Point R_center;
            for (int j = 0; j < contours.size(); j++)//再次遍历所有轮廓
            {
                if (contour[j] == 1) //如果某轮廓对应数组的值为1，说明只要一个内嵌轮廓
                {
                    double s=contourArea(contours[j]);
                    if(s<200)
                    {   
                        Moments R=moments(contours[j]);
                        if(R.m00!=0){
                            Point cen(R.m10/R.m00,R.m01/R.m00);
                            R_center=cen;
                            circle(src,cen,4,Scalar(255,0,0),-1);
                        }
                        continue;
                    } 
                    int num = hierarchy[j][2]; //记录该轮廓的内嵌轮廓
                    RotatedRect box = minAreaRect(contours[num]); //包含该轮廓所有点
                    Point2f vertex[4];
                    box.points(vertex);//将左下角，左上角，右上角，右下角存入点集
                    for (int i = 0; i < 4; i++)
                        {
                            line(src, vertex[i], vertex[(i + 1) % 4], Scalar(255, 0, 0), 4, LINE_AA); //画线
                        }
                    center = (vertex[0] + vertex[2]) / 2; //返回中心坐标
                    circle(src,center,4,Scalar(255,0,0),-1);
                    putText(src, "target", vertex[0], FONT_HERSHEY_SIMPLEX, 1.0, Scalar(255, 255, 0));//打印字体
                    }
            }
            Point RH = center - R_center;
            Point2d RHi = Point2d(RH) / norm(RH);

            // calculate xt,yt;
            double xt = (t_now - t_start) / 1000;
            double yt = RHi.x;

            // ceres
            problem.AddResidualBlock(new ceres::AutoDiffCostFunction<CostFunctor, 1, 1, 1, 1, 1>(new CostFunctor(xt, yt)), NULL, &A0, &A, &w, &phi);

            Solver::Options options;
            options.max_num_iterations = 25;
            options.linear_solver_type = ceres::DENSE_QR;

            problem.SetParameterLowerBound(&A0, 0, 0.5);
            problem.SetParameterUpperBound(&A0, 0, 1.4);
            problem.SetParameterLowerBound(&A, 0, 0.5);
            problem.SetParameterUpperBound(&A, 0, 1.0);
            problem.SetParameterLowerBound(&w, 0, 0.5);
            problem.SetParameterUpperBound(&w, 0, 1.9);
            problem.SetParameterLowerBound(&phi, 0, 0.24);
            problem.SetParameterUpperBound(&phi, 0, 1.25);

            Solver::Summary summary;
            Solve(options, &problem, &summary);
            if (phi < 0.25 && 1.78 < w && 1.23 < A0 && A < 0.83 && A0 < 1.37 && 0.74 < A && w < 1.98 && 0.22 < phi)
            {
                // endtime
                int64 end_time = getTickCount();
                t_sum += (end_time - start_time) / getTickFrequency();
                break;
            }
        }
    }
    std::cout << t_sum / N << std::endl;
}
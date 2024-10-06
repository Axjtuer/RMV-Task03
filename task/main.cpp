#include "windmill.hpp"
#include "ceres/ceres.h"
#include "glog/logging.h"

using namespace std;
using namespace cv;
using namespace ceres;

//计算拟合残差
struct CosRes
{
    CosRes(double t, double y) : t_(t), y_(y) {}

    template <typename T>
    bool operator()(const T *const A, const T *const w, const T *const phi, const T *const b, T *res) const
    {
        T pre_y=ceres::cos(*A / *w * (ceres::sin(*w * t_ + *phi)-ceres::sin(*phi))+*b * t_);
        *res = y_ - pre_y;
        return true;
    }

private:
    const double t_;
    const double y_;
};

int main()
{
    double t_sum = 0;
    const int N = 10;
    for (int num = 0; num < N; num++)
    {
        std::chrono::milliseconds t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
        double t_start = (double)t.count();
        WINDMILL::WindMill wm(t_start);
        Mat src;

        Problem problem;
        double A = 1.785;
        double w = 0.884;
        double phi = 0.24;
        double b = 0.305;

        // starttime
        int64 start_time = getTickCount();

        while (1)
        {
            t = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch());
            double t_now = (double)t.count();
            src = wm.getMat(t_now);

            //灰度化+二值化
            Mat grayimg;
            cvtColor(src,grayimg,COLOR_BGR2GRAY);
            Mat binary;
            threshold(grayimg,binary,75,225,THRESH_BINARY);
            //处理轮廓
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
            
            double time=(t_now-t_start) / 1000;
            double angle=(center.x-R_center.x)/norm(center-R_center);
            //使用ceres Solver拟合
            problem.AddResidualBlock(new AutoDiffCostFunction<CosRes, 1, 1, 1, 1, 1>(new CosRes(time, angle)), nullptr,&A, &w, &phi,&b);

            Solver::Options options;
            options.max_num_iterations=25;
            options.linear_solver_type=DENSE_QR;

            problem.SetParameterLowerBound(&A, 0, 0.4);
            problem.SetParameterUpperBound(&A, 0, 1.1);
            problem.SetParameterLowerBound(&w, 0, 0.45);
            problem.SetParameterUpperBound(&w, 0, 1.89);
            problem.SetParameterLowerBound(&phi, 0, 0.24);
            problem.SetParameterUpperBound(&phi, 0, 1.15);
            problem.SetParameterLowerBound(&b, 0, 0.49);
            problem.SetParameterUpperBound(&b, 0, 1.3);

            Solver::Summary summary;
            Solve(options, &problem, &summary);
            if (0.74575 < A && A < 0.82425 &&
                1.7898 < w && w < 1.9782 &&
                0.228 < phi && phi < 0.252 &&
                1.23975 < b && b < 1.37025)
            {
                long long end_time=getTickCount();
                t_sum+=(end_time-start_time)/getTickFrequency();
                break;
            }
        }
    }
    std::cout << t_sum / N << std::endl;
}
        

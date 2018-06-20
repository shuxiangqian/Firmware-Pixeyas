/*
 * This program aims to estimate the location of indoor UAV under the prior map,
 four mutually perpendicular srf01 sonars and IMU data.
 * @author Xqc
 */

//#include <algorithm>
#include <px4_posix.h>
#include <mathlib/mathlib.h>
#include <uORB/uORB.h>
#include <uORB/topics/control_state.h>
#include <uORB/topics/sonar_distance.h>
#include <drivers/drv_hrt.h>
#include <uORB/topics/ekf_localization.h>
#include <uORB/topics/distance_sensor.h>	// mb12xx
//#include <uORB/topics/sonar_distance_theory.h>

// used for print; uncomment when no longer needed
//double dist[4];
//float RAY[4], WALL[4];
//math::Matrix<2, 1> input;
//math::Vector<3> ATT;
//double H[4][4];
//uint64_t time1, time2, dt;
// --------------print---------------------

static bool task_should_exit = false;
static bool task_running = false;
static int control_task;
float _yaw;
struct control_state_s sensor;
struct distance_sensor_s distance;
int _ctrl_state_sub; 						//control state subscription
int _distance_sensor_sub;          		// four mb12xx sonar sensor data

extern "C" __EXPORT int mc_localization_EKF_main(int argc, char *argv[]);
int task_main(int argc, char *argv[]);

//以下宏定义为室内定位算法需要用到的调节参数
#define PI							3.1415926
#define SONAR_NUMBER_SET 			4				//超声波数量设置

#define RAY_NUMBER_SET   			5				//超声波射线模型中射线数，目前可选参数包括1/3/5/7/21
#define MAP_NUMBER_SET   			13				//先验地图维数，目前可选参数包括4/13/25

struct PRIOR_MAP {
	float x[MAP_NUMBER_SET] = { 0, 5, 5, 9, 9, 14, 14, 12, 12, 7, 7, 0, 0 };
	float y[MAP_NUMBER_SET] = { 0, 0, 2, 2, 0, 0, 4, 4, 10, 10, 5, 5, 0 };
};
PRIOR_MAP _map_test;                          	// the prior map used for test

struct PRIOR_SONAR {
	float x[RAY_NUMBER_SET] = { 1.00, 3.50, 6.00, 1.00, 3.50 };
	float y[RAY_NUMBER_SET] = { -0.2, -0.1375, 0, 0.2, 0.1375 };
};
PRIOR_SONAR _sonar_model;                    	// the srf01 sonar sensor model

float Wall_Angle[MAP_NUMBER_SET - 1] = { 0, PI / 2, 0, PI / 2, 0, PI / 2, 0, PI
		/ 2, 0, PI / 2, 0, PI / 2 };
//
struct POINTF {
	POINTF(float a, float b) :
			x(a), y(b) {
	}
	POINTF() {
	}

	float x;
	float y;
};
POINTF _map_start_point, // the map start point for calculate the intersection point
		_map_end_point, // the map end point for calculate the intersection point
		_sonar_start_point,                   // the sonar model start point
		_sonar_end_point,                       // the sonar model end point
		_sonar_end_point_rotation, // the sonar model end point after rotation
		_intersection_point; // the intersection point between the sonar model and the prior map

// 使用提示函数
static void usage(const char *reason);
int sonar_data_update(bool force); // aim to find sonar data status,if update,we accept it
int sensor_data_update(bool force); // aim to find sensor data status,if update,we accept it

void calculateline(POINTF p1, POINTF p2, float &a, float &b, float &c); //	calculate the line on the point p1 and p2
bool Equal(float f1, float f2);           // judge f1 and f2 is equal or not
bool Bigger(const POINTF &p1, const POINTF &p2); // judge the point p1 is bigger than p2 or not
float Cross_product(const POINTF &p1, const POINTF &p2);          // 计算两向量外积
void Swap(float &f1, float &f2);              // swap the value of f1 and f2
void Swap_struct(POINTF &p1, POINTF &p2);  // swap the position of p1 and p2

void ctrlStateUpdate(bool&);
void distanceSensorUpdate(bool&);

//判定两线段位置关系，并求出交点(如果存在)
float Intersection(POINTF p1, POINTF p2, POINTF p3, POINTF p4,
		POINTF &intersection_point);

//void Intersection(POINTF, POINTF, double, POINTF, POINTF, POINTF &, float&);

// the purpose of this function is to get the the minimum distance between the cross point and the sonar model origin
void sonar_value_minimum(POINTF location, float yaw,
		math::Matrix<4, 10> &sonar_map_label,
		math::Matrix<4, 1> &distance_theory_minimum, float[], float[]);

//计算雅克比矩阵
//此雅克比矩阵通过几何方法计算，目的是节省计算量
void caculateJacobi(math::Matrix<4, 1> &, float[], float[],
		math::Matrix<4, 4> &);

//此函数用于计算直线表达式，已知两点p1和p2，得到直线ax+by+c
void calculateline(POINTF p1, POINTF p2, float &a, float &b, float &c) {
	a = p2.y - p1.y;
	b = p1.x - p2.x;
	c = p2.x * p1.y - p1.x * p2.y;
}

// judge f1 and f2 is equal or not
bool Equal(float f1, float f2) {
	return (fabs(f1 - f2) < 1e-4);
}

//比较两点坐标大小，先比较x坐标，若相同则比较y坐标
bool Bigger(const POINTF &p1, const POINTF &p2) {
	return (p1.x > p2.x || (Equal(p1.x, p2.x) && p1.y > p2.y));
}

//计算两向量外积
float Cross_product(const POINTF &p1, const POINTF &p2) {
	return (p1.x * p2.y - p1.y * p2.x);
}
//swap f1 and f2
void Swap(float &f1, float &f2) {
	float temp;
	temp = f1;
	f1 = f2;
	f2 = temp;
}
//swap the position p1 and p2
void Swap_struct(POINTF &p1, POINTF &p2) {
	Swap(p1.x, p2.x);
	Swap(p1.y, p2.y);
}

//void Intersection(POINTF _sonar_start, POINTF _sonar_end, double angle,
//		POINTF _map_start, POINTF _map_end, POINTF &intersection_point,
//		float& dist1) {
//
//	POINTF sonar1, v2(0, 0), v1(0, 0), v3;
//	float distr_ray;
//	//求方向向量
//	sonar1 = POINTF(cos(angle), sin(angle));
//	v1 = POINTF(_map_start.x - _sonar_start.x, _map_start.y - _sonar_start.y);
//	v2 = POINTF(_map_end.x - _sonar_start.x, _map_end.y - _sonar_start.y);
//	//计算lambda
//	float lamb = -(v1.y * sonar1.x + v1.x * sonar1.y)
//			/ (v2.x * sonar1.y + v2.y * sonar1.x);
////    cout << "lamb = " << lamb << endl;
//	//计算交点
//	if (lamb >= 0) {
//		intersection_point.x = (_map_start.x + lamb * _map_end.x) / (1 + lamb);
//		intersection_point.y = (_map_start.y + lamb * _map_end.y) / (1 + lamb);
//		v3 = POINTF(intersection_point.x - _sonar_start.x,
//				intersection_point.y - _sonar_start.y);
//		dist1 = sqrtf(v3.x * v3.x + v3.y * v3.y);
//
//		distr_ray = sqrtf(
//				(_sonar_end.x - _sonar_start.x)
//						* (_sonar_end.x - _sonar_start.x)
//						+ (_sonar_end.y - _sonar_start.y)
//								* (_sonar_end.y - _sonar_start.y));
//
//		if (dist1 > distr_ray + 0.0001f) {
//			intersection_point.x = 1000;
//			intersection_point.y = 1000;
//		}
//
//		if (sonar1.x * v3.x + sonar1.y * v3.y < 0) {
//			intersection_point.x = 1000;
//			intersection_point.y = 1000;
//		}
//
//	} else {
//		//no intersection.
//	}
//}
//判定两线段位置关系，并求出交点(如果存在)。返回值交于线上(2)，正交(1)，无交(0)
float Intersection(POINTF p1, POINTF p2, POINTF p3, POINTF p4,
		POINTF &intersection_point) {
	// 为方便运算，保证各线段的起点在前，终点在后。
	if (Bigger(p1, p2)) {
		Swap_struct(p1, p2);
	}
	if (Bigger(p3, p4)) {
		Swap_struct(p3, p4);
	}
	// 将线段按起点坐标排序。若线段1的起点较大，则将两线段交换
	if (Bigger(p1, p3)) {
		Swap_struct(p1, p3);
		Swap_struct(p2, p4);
	}
	//计算向量及其外积
	POINTF v1 = { p2.x - p1.x, p2.y - p1.y }, v2 = { p4.x - p3.x, p4.y - p3.y };
	float Corss = Cross_product(v1, v2);
	//先进行快速排斥试验
	//x坐标已有序，可直接比较。y坐标要先求两线段的最大和最小值
	float ymax1 = p1.y, ymin1 = p2.y, ymax2 = p3.y, ymin2 = p4.y;
	if (ymax1 < ymin1) {
		Swap(ymax1, ymin1);
	}
	if (ymax2 < ymin2) {
		Swap(ymax2, ymin2);
	}
	//如果以两线段为对角线的矩形不相交，则无交点
	if (p1.x > p4.x || p2.x < p3.x || ymax1 < ymin2 || ymin1 > ymax2) {
		return 0;
	}
	//下面进行跨立试验
	POINTF vs1 = { p1.x - p3.x, p1.y - p3.y }, vs2 =
			{ p2.x - p3.x, p2.y - p3.y };
	POINTF vt1 = { p3.x - p1.x, p3.y - p1.y }, vt2 =
			{ p4.x - p1.x, p4.y - p1.y };
	float s1v2, s2v2, t1v1, t2v1;
	//根据外积结果判定否交于线上
	s1v2 = Cross_product(vs1, v2);
	s2v2 = Cross_product(vs2, v2);
	t1v1 = Cross_product(vt1, v1);
	t2v1 = Cross_product(vt2, v1);
	//未交于线上，则判定是否相交
	if (s1v2 * s2v2 > 0 || t1v1 * t2v1 > 0) {
		return 0;
	}
	//以下为相交的情况,计算二阶行列式的两个常数项
	float ConA = p1.x * v1.y - p1.y * v1.x;
	float ConB = p3.x * v2.y - p3.y * v2.x;
	//计算行列式D1和D2的值，除以系数行列式的值，得到交点坐标
	intersection_point.x = (ConB * v1.x - ConA * v2.x) / Corss;
	intersection_point.y = (ConB * v1.y - ConA * v2.y) / Corss;
	//正交返回1
	return 1;
}

void sonar_value_minimum(POINTF location, float yaw,
		math::Matrix<4, 10> &sonar_map_label,
		math::Matrix<4, 1> &distance_theory_minimum, float ray_angle[],
		float wall_angle[]) {

	float _distance_theory[4];
	float para_a;
	float para_b;
	float para_c;
	float side_1;
	float side_2;
	float side;

	_distance_theory[0] = 7;
	_distance_theory[1] = 7;
	_distance_theory[2] = 7;
	_distance_theory[3] = 7;


	ray_angle[0] = -1;
	ray_angle[1] = -1;
	ray_angle[2] = -1;
	ray_angle[3] = -1;

	wall_angle[0] = -1;
	wall_angle[1] = -1;
	wall_angle[2] = -1;
	wall_angle[3] = -1;

	/*
	 sonar_map_label用于存储有交点的超声波射线和地图线段，首先全部赋值-1
	 当对应位置有交点时，则以交点编号和位置数据替换
	 如果没有交点，则为-1不变，作为后续判断的标志位
	 */
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 10; j++) {
			sonar_map_label(i, j) = -1;
		}
	}

	float _distance_theory_temp = 0;
	float cos_value;
	float sin_value;
	float sin_yaw;
	float cos_yaw;
	sin_yaw = sin(yaw);
	cos_yaw = cos(yaw);

	_sonar_start_point.x = location.x; //Initialize the sonar model position
	_sonar_start_point.y = location.y; //Initialize the sonar model position

	for (int i = 0; i < SONAR_NUMBER_SET; i++) {
		cos_value = cos(90 * i * PI / 180);
		sin_value = sin(90 * i * PI / 180);

		for (int j = 0; j < RAY_NUMBER_SET; j++) {
			_sonar_end_point_rotation.x = (_sonar_model.x[j] * cos_yaw)
					- (_sonar_model.y[j] * sin_yaw);
			_sonar_end_point_rotation.y = (_sonar_model.x[j] * sin_yaw)
					+ (_sonar_model.y[j] * cos_yaw);

			//the ith sonar model position x'=x*cos(90*(i-1))-y*sin(90*(i-1))
			_sonar_end_point.x = location.x
					+ (_sonar_end_point_rotation.x * cos_value)
					- (_sonar_end_point_rotation.y * sin_value);
			//the ith sonar model position y'=x*sin(90*(i-1))+y*cos(90*(i-1))
			_sonar_end_point.y = location.y
					+ (_sonar_end_point_rotation.x * sin_value)
					+ (_sonar_end_point_rotation.y * cos_value);


//			angle = atan2((_sonar_end_point.y - _sonar_start_point.y),
//					(_sonar_end_point.x - _sonar_start_point.x));
			//由点_sonar_start_point和点_sonar_end_point计算得到穿过两点的直线para_a*x + para_b*y + para_c
			calculateline(_sonar_start_point, _sonar_end_point, para_a, para_b,
					para_c);

			_distance_theory_temp = 10;

			for (int k = 1; k < MAP_NUMBER_SET; k++) {


				_map_start_point.x = _map_test.x[k - 1]; //the start point of map line
				_map_start_point.y = _map_test.y[k - 1]; //the start point of map line
				_map_end_point.x = _map_test.x[k]; //the end point of map line
				_map_end_point.y = _map_test.y[k]; //the end point of map line

				_intersection_point.x = 1000;
				_intersection_point.y = 1000;

				//判断点_map_start_point和点_map_end_point与直线para_a*x + para_b*y + para_c的位置关系
				//当两点在直线同一侧时，side>0；当两点在直线两侧时，side<0；当至少一个点位于直线上时，side=0.
				side_1 = para_a * _map_start_point.x
						+ para_b * _map_start_point.y + para_c;
				side_2 = para_a * _map_end_point.x + para_b * _map_end_point.y
						+ para_c;
				side = side_1 * side_2;

				if (side > 0) {

				} else {
					// for calculate the intersection between the sonar model and the prior map
					Intersection(_sonar_start_point, _sonar_end_point,
							_map_start_point, _map_end_point,
							_intersection_point);

//					Intersection(_sonar_start_point, _sonar_end_point, angle,
//							_map_start_point, _map_end_point,
//							_intersection_point, _distance_theory_temp);
				}
				// for calculate the distance between the intersection point and the estimation point(the sonar model origin)
				if (((int) _intersection_point.x == 1000)
						&& ((int) _intersection_point.y == 1000)) {

				} else {
					_distance_theory_temp = sqrtf(
							(_intersection_point.x - location.x)
									* (_intersection_point.x - location.x)
									+ (_intersection_point.y - location.y)
											* (_intersection_point.y
													- location.y));

					if (_distance_theory[i] > _distance_theory_temp) {
						//当有交点时开始赋值相应数据
						//i表示第i个超声波
						sonar_map_label(i, 0) = j;					//第j根射线
						sonar_map_label(i, 1) = _sonar_start_point.x;//第j根射线起始位置X坐标
						sonar_map_label(i, 2) = _sonar_start_point.y;//第j根射线起始位置Y坐标
						sonar_map_label(i, 3) = _sonar_end_point.x;	//第j根射线结束位置X坐标
						sonar_map_label(i, 4) = _sonar_end_point.y;	//第j根射线结束位置Y坐标
						sonar_map_label(i, 5) = k;					//第k条地图线段
						sonar_map_label(i, 6) = _map_start_point.x;	//第k条地图线段起始位置X坐标
						sonar_map_label(i, 7) = _map_start_point.y;	//第k条地图线段起始位置Y坐标
						sonar_map_label(i, 8) = _map_end_point.x;//第k条地图线段结束位置X坐标
						sonar_map_label(i, 9) = _map_end_point.y;//第k条地图线段结束位置Y坐标

						_distance_theory[i] = _distance_theory_temp;

						float dy = _sonar_end_point.y - _sonar_start_point.y;
						float dx = _sonar_end_point.x - _sonar_start_point.x;

						ray_angle[i] = atan2(dy, dx) * (180 / PI);
						if (dx < 0 && dy < 0) {
							ray_angle[i] += 360;
						} else if (dx > 0 && dy < 0) {
							ray_angle[i] += 360;
						}
						wall_angle[i] = Wall_Angle[k - 1] * (float) (180 / PI);
//						break;
					}
				}
			}
			usleep(2);
		}
	}
	// store the minimum distance in distance_theory_minimum
	distance_theory_minimum(0, 0) = _distance_theory[0];
	distance_theory_minimum(1, 0) = _distance_theory[1];
	distance_theory_minimum(2, 0) = _distance_theory[2];
	distance_theory_minimum(3, 0) = _distance_theory[3];
}

void caculateJacobi(math::Matrix<4, 1> &sonarTheory, float rayAngle[],
		float wallAngle[], math::Matrix<4, 4> &h) {

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			h(i, j) = 0;

	for (int i = 0; i < 4; i++) {
		if (fabs(sonarTheory(i, 0) - 7) > 0.001) {
			h(i, 0) = -sin(wallAngle[i] * (float) (PI / 180))
					/ sin((wallAngle[i] - rayAngle[i]) * (float) (PI / 180));
			h(i, 2) = cos(wallAngle[i] * (float) (PI / 180))
					/ sin((wallAngle[i] - rayAngle[i]) * (float) (PI / 180));
		}
	}
}

void ctrlStateUpdate(bool& updated) {
	math::Matrix<3, 3> _R; 			// rotation matrix from attitude quaternions
	orb_check(_ctrl_state_sub, &updated);
	if (updated) {
		orb_copy(ORB_ID(control_state), _ctrl_state_sub, &sensor);
		math::Quaternion q_att(sensor.q[0], sensor.q[1], sensor.q[2],
				sensor.q[3]);
		_R = q_att.to_dcm();
		math::Vector<3> euler_angles;
		euler_angles = _R.to_euler();
		_yaw = euler_angles(2);		// the 3rd element of eulerian angle is yaw
//		ATT = euler_angles;
//		ATT = ATT*(180/PI);
	}
}

void distanceSensorUpdate(bool& updated) {
	orb_check(_distance_sensor_sub, &updated);
	if (updated) {
		orb_copy(ORB_ID(distance_sensor), _distance_sensor_sub, &distance);
	}
}
//任务主函数
int task_main(int argc, char *argv[]) {
	usleep(1000);
	warnx("mc_localization_EKF is successful!\n");
	//线程启动标志
	task_running = true;

	warnx("mc_localization_EKF start successfully\n");

// Load the prior map，加载地图模型

	if (MAP_NUMBER_SET == 4) {
		//4维的地图
		auto map_data_boundary_x = std::initializer_list<float>(
				{ 0.00, 6, 6, 0 });
		std::copy(map_data_boundary_x.begin(), map_data_boundary_x.end(),
				_map_test.x);
		auto map_data_boundary_y = std::initializer_list<float>( { 0.00, 0.00,
				7.5, 7.5 });
		std::copy(map_data_boundary_y.begin(), map_data_boundary_y.end(),
				_map_test.y);
	}

	if (MAP_NUMBER_SET == 13) {
		//13维的地图
		auto map_data_boundary_x = std::initializer_list<float>( { 0, 6.30,
				6.30, 0, 0, -1.95, -1.95, -11, -11, -1.95, -1.95, 0, 0 });
		std::copy(map_data_boundary_x.begin(), map_data_boundary_x.end(),
				_map_test.x);

		auto map_data_boundary_y = std::initializer_list<float>( { 0, 0, 7.26,
				7.26, 14.26, 14.26, 5.23, 5.23, 2.03, 2.03, -7.00, -7.00, 0 });
		std::copy(map_data_boundary_y.begin(), map_data_boundary_y.end(),
				_map_test.y);
	}

	if (MAP_NUMBER_SET == 25) {
		//25维的地图
		auto map_data_boundary_x = std::initializer_list<float>( { 0, 6.30,
				6.30, 0, 0, -1.95, -1.95, -11, -11, -9.75, -9.75, -8.60, -8.60,
				-7.10, -7.10, -5.95, -5.95, -4.45, -4.45, -3.30, -3.30, -1.95,
				-1.95, 0, 0 });
		std::copy(map_data_boundary_x.begin(), map_data_boundary_x.end(),
				_map_test.x);

		auto map_data_boundary_y = std::initializer_list<float>( { 0, 0, 7.26,
				7.26, 14.26, 14.26, 5.23, 5.23, 2.03, 2.03, 1.48, 1.48, 2.03,
				2.03, 1.48, 1.48, 2.03, 2.03, 1.48, 1.48, 2.03, 2.03, -7.00,
				-7.00, 0 });
		std::copy(map_data_boundary_y.begin(), map_data_boundary_y.end(),
				_map_test.y);
	}

// Load the srf01 sonar model，加载超声波射线模型

	if (RAY_NUMBER_SET == 1) {
		//1根射线模型
		auto sonar_data_x = std::initializer_list<float>( { 6.00 });
		std::copy(sonar_data_x.begin(), sonar_data_x.end(), _sonar_model.x);
		auto sonar_data_y = std::initializer_list<float>( { 0.00 });
		std::copy(sonar_data_y.begin(), sonar_data_y.end(), _sonar_model.y);
	}

	if (RAY_NUMBER_SET == 3) {
		//3根射线模型
		auto sonar_data_x = std::initializer_list<float>( { 1.00, 6.00, 1.00 });
		std::copy(sonar_data_x.begin(), sonar_data_x.end(), _sonar_model.x);
		auto sonar_data_y = std::initializer_list<float>( { -0.2, 0.00, 0.2 });
		std::copy(sonar_data_y.begin(), sonar_data_y.end(), _sonar_model.y);
	}

	if (RAY_NUMBER_SET == 5) {
		//5根射线模型
		auto sonar_data_x = std::initializer_list<float>( { 1.00, 3.50, 6.00,
				1.00, 3.50 });
		std::copy(sonar_data_x.begin(), sonar_data_x.end(), _sonar_model.x);
		auto sonar_data_y = std::initializer_list<float>( { -0.2, -0.1375, 0,
				0.2, 0.1375 });
		std::copy(sonar_data_y.begin(), sonar_data_y.end(), _sonar_model.y);
	}

	if (RAY_NUMBER_SET == 7) {
		//7根射线模型
		auto sonar_data_x = std::initializer_list<float>( { 1.00, 2.00, 3.50,
				6.00, 1.00, 2.00, 3.50 });
		std::copy(sonar_data_x.begin(), sonar_data_x.end(), _sonar_model.x);
		auto sonar_data_y = std::initializer_list<float>( { -0.2, -0.135,
				-0.1375, 0, 0.2, 0.135, 0.1375 });
		std::copy(sonar_data_y.begin(), sonar_data_y.end(), _sonar_model.y);
	}

	if (RAY_NUMBER_SET == 21) {
		//21根射线模型
		auto sonar_data_x = std::initializer_list<float>( { 0.60, 1.00, 1.50,
				2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 0.60, 1.00,
				1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00 });
		std::copy(sonar_data_x.begin(), sonar_data_x.end(), _sonar_model.x);

		auto sonar_data_y = std::initializer_list<float>( { -0.015, -0.2, -0.19,
				-0.135, -0.15, -0.1575, -0.1375, -0.03, -0.04, -0.005, 0, 0.015,
				0.2, 0.19, 0.135, 0.15, 0.1575, 0.1375, 0.03, 0.04, 0.005 });
		std::copy(sonar_data_y.begin(), sonar_data_y.end(), _sonar_model.y);
	}
	

	_ctrl_state_sub = orb_subscribe(ORB_ID(control_state));  // subscribe the control state
	_sonar_sub_fd = orb_subscribe(ORB_ID(sonar_distance));   // subscribe sonar data
	// Initialize the estimation location
	math::Matrix<2, 1> point_estimation;
	point_estimation(0, 0) = POINT_ESTIMATION_X_SET;
	point_estimation(1, 0) = POINT_ESTIMATION_Y_SET;
	math::Matrix<2, 1> point_previous;
	point_previous(0, 0) = POINT_PREVIOUS_X_SET;
	point_previous(1, 0) = POINT_PREVIOUS_Y_SET;
	math::Matrix<4, 1> sonarvalue_reality;
	math::Matrix<4, 10> sonar_map_label;
	//任务循环
	while (!task_should_exit) {
		usleep(5000);
		
	warnx("mc_localization_EKF start successfully\n");

	//rc_channels
	//orb_subscribe(ORB_ID(rc_channels));
	_ctrl_state_sub = orb_subscribe(ORB_ID(control_state)); // subscribe the control state
//	_distance_sensor_sub = orb_subscribe(ORB_ID(distance_sensor)); // subscribe sonar data
	struct ekf_localization_s mav_position;
	orb_advert_t ekf_localization_pub = orb_advertise(ORB_ID(ekf_localization),
			&mav_position);

	// Initialize the estimation location
	float sonar_ray_angle[4];
	float sonar_wall_angle[4];
	float T = 0.008f;		// imu sample
	math::Matrix<4, 4> Jacobi;
	math::Matrix<4, 10> sonar_map_label;

	// position corrected
	math::Matrix<4, 1> X;
	X(0, 0) = 1.0f;
	X(1, 0) = 0.0f;
	X(2, 0) = 1.0f;
	X(3, 0) = 0.0f;
	// x_hat = A*x+B*u
	math::Matrix<4, 1> X_hat;
	X_hat.zero();
	// state input
	math::Matrix<2, 1> u;
	math::Matrix<4, 4> A;
	A.zero();
	A(0, 0) = 1;
	A(0, 1) = T;
	A(1, 1) = 1;
	A(2, 2) = 1;
	A(2, 3) = T;
	A(3, 3) = 1;

	math::Matrix<4, 2> B;
	B(0, 0) = T * T / 2;
	B(0, 1) = 0;
	B(1, 0) = T;
	B(1, 1) = 0;
	B(2, 0) = 0;
	B(2, 1) = T * T / 2;
	B(3, 0) = 0;
	B(3, 1) = T;

	math::Matrix<4, 4> P;		//var
	P.identity();
	P = P * 200;
	math::Matrix<4, 4> P0;

	math::Matrix<4, 4> Q;	// process noise
	Q.zero();
	Q(0, 0) = 1;
	Q(1, 1) = 0.2f;
	Q(2, 2) = 1;
	Q(3, 3) = 0.2f;

	math::Matrix<4, 4> R;	// process noise
	R.identity();
	R = R * 0.007f;

	math::Matrix<4, 4> I;
	I.identity();
//
	math::Matrix<4, 4> K;
	K.zero();

	math::Matrix<4, 1> sonar_measure;

	math::Matrix<4, 1> sonarvalue_theory;

	//任务循环
	while (!task_should_exit) {

		bool updated = true;
		ctrlStateUpdate(updated);
		if (updated) {
			u(0, 0) = sensor.x_acc * (float) cos(_yaw)
					+ sensor.y_acc * (float) sin(_yaw);
			u(1, 0) = sensor.x_acc * (float) sin(_yaw)
					- sensor.y_acc * (float) cos(_yaw);
//			input(0, 0) = sensor.x_acc;
//			input(1, 0) = sensor.y_acc;
			// 1
			X_hat = A * X + B * u;
//			// 2
			P0 = A * P * A.transposed() + Q;
			distanceSensorUpdate(updated);
			if (updated) {
				for (int i = 0; i < 4; i++)
					sonar_measure(i, 0) = distance.distance[i];

				POINTF position_estimation(X_hat(0, 0), X_hat(2, 0));
//			POINTF position_estimation = { 4, 3 };
//			_yaw = 0 * (PI / 180);

				sonar_value_minimum(position_estimation, _yaw, sonar_map_label,
						sonarvalue_theory, sonar_ray_angle, sonar_wall_angle);

				caculateJacobi(sonarvalue_theory, sonar_ray_angle,
						sonar_wall_angle, Jacobi);

				K = P0 * Jacobi.transposed()
						* ((Jacobi * P0 * Jacobi.transposed() + R).inversed());
				X = X_hat + K * (sonar_measure - sonarvalue_theory);
				P = (I - K * Jacobi) * P0;
			} else {
				X = X_hat;
				P = P0;
			}
			mav_position.x = X(0, 0);
			mav_position.vx = X(1, 0);
			mav_position.y = X(2, 0);
			mav_position.vy = X(3, 0);

//			mav_position.x = 1;
//			mav_position.vx = 2;
//			mav_position.y = 3;
//			mav_position.vy = 4;

			orb_publish(ORB_ID(ekf_localization), ekf_localization_pub,
					&mav_position);
		}

		POINTF position_estimation = { 4, 3 };
//		time1 = hrt_absolute_time() / 1000;
		sonar_value_minimum(position_estimation, _yaw, sonar_map_label,
				sonarvalue_theory, sonar_ray_angle, sonar_wall_angle);
		caculateJacobi(sonarvalue_theory, sonar_ray_angle, sonar_wall_angle,
				Jacobi);
//		time2 = hrt_absolute_time() / 1000;
//		dt = time2 - time1;

//		 ----------------print---------------
//		dist[0] = sonarvalue_theory(0, 0);
//		dist[1] = sonarvalue_theory(1, 0);
//		dist[2] = sonarvalue_theory(2, 0);
//		dist[3] = sonarvalue_theory(3, 0);
//		RAY[0] = sonar_ray_angle[0];
//		RAY[1] = sonar_ray_angle[1];
//		RAY[2] = sonar_ray_angle[2];
//		RAY[3] = sonar_ray_angle[3];
//		WALL[0] = sonar_wall_angle[0];
//		WALL[1] = sonar_wall_angle[1];
//		WALL[2] = sonar_wall_angle[2];
//		WALL[3] = sonar_wall_angle[3];

//		for(int i = 0;i<4;i++)
//			for(int j=0;j<4;j++)
//				H[i][j] = Jacobi(i, j);
		//-----------------print---------------
		usleep(10000);
	}
	//关闭线程，串口
	warnx("pos_estimator_sonar_imu exiting.\n");
	task_running = false;

	return 0;
}

int mc_localization_EKF_main(int argc, char *argv[]) {
	if (argc < 2) {
		usage("[pos_estimator_sonar_imu]missing command");
		return 1;
	}
	if (!strcmp(argv[1], "start")) {
		if (task_running) {
			warnx("[pos_estimator_sonar_imu]already running\n");
			return 0;
		}
		task_should_exit = false;

		control_task = px4_task_spawn_cmd("mc_localization_EKF", SCHED_DEFAULT,
		SCHED_PRIORITY_MAX - 5, 4000, task_main,

				(argv) ? (char * const *) &argv[2] : (char * const *) NULL);
		return (0);
	}
	if (!strcmp(argv[1], "stop")) {
		task_should_exit = true;
		return (0);
	}
	if (!strcmp(argv[1], "status")) {
		if (task_running) {
			warnx("[mc_localization_EKF] running");
			while (1) {
//				printf(
//						"[EKF] dist[0] = %.2f, dist[1] = %.2f, dist[2] = %.2f, dist[3] = %.2f \n",
//						dist[0], dist[1], dist[2], dist[3]);
//				printf("[EKF] time = %lld \n", dt);
//				printf("[EKF] ray_angle = %.2f, %.2f, %.2f %.2f \n",
//						(double) RAY[0], (double) RAY[1], (double) RAY[2],
//						(double) RAY[3]);
//				printf("[EKF] wall_angle = %.2f, %.2f, %.2f %.2f \n",
//						(double) WALL[0], (double) WALL[1], (double) WALL[2],
//						(double) WALL[3]);
//				printf("input = \n");
//				input.print();
//				printf("\n");
//				printf("attitude = \n");
//				ATT.print();
//				printf("\n");
//				printf("[EKF] H = \n");
//				for (int i = 0; i < 4; i++) {
//					for (int j = 0; j < 4; j++) {
//						printf("%.2f  ", H[i][j]);
//					}
//					printf("\n");
//				}
				printf("\n");
				printf("==============press CTRL+C to abort==============\n");
				char c;
				struct pollfd fds;
				int ret;
				fds.fd = 0;
				fds.events = POLLIN;
				ret = poll(&fds, 1, 0);
				if (ret > 0) {
					read(0, &c, 1);
					if (c == 0x03 || c == 0x63 || c == 'q') {
						warnx("User abort\n");
						break;
					}
				}
				usleep(400000);		//500ms
			}
		} else {
			warnx("[pos_estimator_sonar_imu]stopped");
		}
		return (0);
	}
	usage("unrecognized command");
	return (1);
}

static void usage(const char *reason) {
	if (reason) {
		fprintf(stderr, "%s\n", reason);
	}
	fprintf(stderr,
			"usage: pos_estimator_sonar_imu {start|stop|status} [param]\n\n");
	exit(1);
}



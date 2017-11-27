/* mc_alt_estimate
 * by xiangqian
 * 2017.8.3     */
//#include <px4_config.h>
//#include <px4_tasks.h>
#include <px4_posix.h>
//#include <algorithm>
//#include <unistd.h>
//#include <stdio.h>
//#include <string.h>
//#include <math.h>
#include <mathlib/mathlib.h>
//#include <arch/board/board.h>
#include <uORB/topics/sonar_distance.h>
#include <uORB/topics/sensor_combined.h>
#include <uORB/topics/alt_estimate.h>
//#include <drivers/drv_hrt.h>
//#include <lib/mathlib/mathlib.h>

extern "C" __EXPORT int mc_alt_estimator_main(int argc, char *argv[]);

//

static bool thread_should_exit = false;
static bool thread_running = false;
static int daemon_task;
// for debug
//uint64_t count12 = 0;
//float acc = 0.0f;
float test[10]={0};
//uint64_t time1,time2,time_cha;
//unsigned int count1=0;

__EXPORT int mc_alt_estimator_main(int argc, char *argv[]);
int mc_alt_estimator_thread_main(int argc, char *argv[]);

int 	_sonar_sub;
int  	_sensor_sub;
// 使用提示函数
static void usage(const char *reason);

static void usage(const char *reason)
{
    if (reason) {
        fprintf(stderr, "%s\n", reason);
    }

    fprintf(stderr, "usage: alt_estimate {start|stop|status} [param]\n\n");
    exit(1);
}

int mc_alt_estimator_main(int argc, char *argv[])
{
    if (argc < 2) {
        usage("[Alt_est] missing command");
        return 1;
    }

    if (!strcmp(argv[1], "start")) {
        if (thread_running) {
            warnx("[Alt_est] already running!\n");
            return 0;
        }

        thread_should_exit = false;
        daemon_task = px4_task_spawn_cmd("mc_alt_estimator",
                                         SCHED_DEFAULT,
                                         SCHED_PRIORITY_MAX - 5,
                                         2000,
										 mc_alt_estimator_thread_main,
                                         (argv) ? (char * const *)&argv[2] : (char * const *)NULL);
        return(0);
    }

    if (!strcmp(argv[1], "stop")) {
        thread_should_exit = true;
        return(0);
    }

    if (!strcmp(argv[1], "status")) {
        if (thread_running) {
            warnx("[Alt_est] running");
            while(1)
            {
            	//printf("count = %lld \n",count12);
            	//warnx("Debug:time=%lld count=%d \n",time_cha,count1);
            	printf("[Test] [0]=%.3f [1]=%.3f [2]=%.3f\n",(double)test[0],(double)test[1],(double)test[2]);
            	printf("[Test] [3]=%.3f [4]=%.3f [5]=%.3f\n",(double)test[3],(double)test[4],(double)test[5]);
            	printf("[Test] [6]=%.3f [7]=%.3f [8]=%.3f\n",(double)test[6],(double)test[7],(double)test[8]);
            	printf("=========Press CTRL+C to abort=========\n");

            	char c;
				struct pollfd fds1;
				int ret;
				fds1.fd=0;
				fds1.events=POLLIN;
				ret=poll(&fds1,1,0);
				if(ret>0)
				{
					read(0,&c,1);
					if(c==0x03||c==0x63||c=='q')
					{
						warnx("User abort\n");
						break;
					}
				}
				usleep(400000);
            }
            return 0;

        } else {
            warnx("[Alt_est] stopped");
        }

        return(0);
    }

    usage("unrecognized command");
    return(1);
}

int mc_alt_estimator_thread_main(int argc, char *argv[])
{
    //线程启动标志
    thread_running = true;
    warnx("[Alt_est] service start successfully!\n");
    _sonar_sub = orb_subscribe(ORB_ID(sonar_distance));
    _sensor_sub = orb_subscribe(ORB_ID(sensor_combined));

    struct sensor_combined_s 	sensor;
    struct sonar_distance_s 	sonar;

    float bias_limit = 0.5f;
    float acc_z = 0.0f, acc_z_pre = 0.0f, T = 0.008f;
    math::Matrix<3,3> Q;	// process noise
    Q(0,0) = 0.05f;  Q(0,1) = 0;  Q(0,2) = 0;
	Q(1,0) = 0;      Q(1,1) = 5;  Q(1,2) = 0;
	Q(2,0) = 0;      Q(2,1) = 0;  Q(2,2) = 0.1f;

	math::Matrix<1,1> R;	// measure noise
	R(0,0) = 0.1f;

	math::Matrix<1,1> m;

    math::Matrix<3,1>  K;		// KF gain
    K(0,0) = 0.0f;
    K(1,0) = 0.0f;
    K(2,0) = 0.0f;

    math::Matrix<3,1> x_hat;
    x_hat(0,0) = 0.0f;
    x_hat(1,0) = 0.0f;
    x_hat(2,0) = -9.8f;

    math::Matrix<3,1> x_hat0;
    x_hat0(0,0) = 0.0f;
	x_hat0(1,0) = 0.0f;
	x_hat0(2,0) = 0.0f;

    math::Matrix<3,3> P;		//var
    P(0,0) = 200; P(0,1) = 0; P(0,2) = 0;
    P(1,1) = 200; P(1,0) = 0; P(1,2) = 0;
    P(2,2) = 200; P(2,0) = 0; P(2,1) = 0;

    math::Matrix<3,3> I;
    I.identity();
    math::Matrix<3,3> P0;

    math::Matrix<3,3> A;
    A(0,0) = 1;  A(0,1) = T;  A(0,2) = 0;
    A(1,0) = 0;  A(1,1) = 1;  A(1,2) = T;
    A(2,0) = 0;  A(2,1) = 0;  A(2,2) = 1;

    math::Matrix<3,3>  Aa;	// A transpose
    Aa = A.transposed();

    math::Matrix<1,3> H;
    H(0,0) = 1; H(0,1) = 0; H(0,2) = 0;

    math::Matrix<3,1> H1;
    H1 = H.transposed();

    math::Matrix<1,1> inv;
    //float inv;

    math::Matrix<3,1> B;
    B(0,0) = 0.0f;  B(1,0) = T;  B(2,0) = 0.0f;

    math::Matrix<1,1> u;	// control input
    math::Matrix<1,1> y;	// measure

    struct alt_estimate_s alt;
    orb_advert_t alt_estimate_pub = orb_advertise(ORB_ID(alt_estimate), &alt);
    //static char count = 0;
    while(!thread_should_exit)
    {
    	bool updated;

    	orb_check(_sensor_sub, &updated);
    	if(updated)			// acc get new data
    	{
    		// test[0] = acc_z_pre;
    		orb_copy(ORB_ID(sensor_combined), _sensor_sub, &sensor);
    		acc_z = sensor.accelerometer_m_s2[2];
    		// 1
    		u(0,0) = -(acc_z_pre+9.88f);
    		x_hat0 = A*x_hat+B*u;
    		// 2
    		P0 = A*P*Aa+Q;

    		orb_check(_sonar_sub, &updated);
    		if(updated)		// new range data
    		{
    			orb_copy(ORB_ID(sonar_distance), _sonar_sub, &sonar);
    			y(0,0) = sonar.distance_filter*0.01f;
    			inv = H*P0*H1;		// 1x3  3x3  3x1 = 1x1
    			inv(0,0) = inv(0,0)+R(0,0);
    			K = P0*H1*inv.inversed();
    			m = H*x_hat0;
    			x_hat = x_hat0+K*(y(0,0)-m(0,0));
    			P = (I-K*H)*P0;
    			alt.flag = true;
    			alt.alt_with_sonar = x_hat(0,0);
    			usleep(1500);
//    			time1=hrt_absolute_time()/1000;
//				time_cha=time1-time2;
//				time2=time1;
//				count1++;
    		}
    		else
    		{
    			x_hat = x_hat0;
    			P = P0;
    			alt.flag = false;
    			usleep(1500);
    		}
    		if(x_hat(2,0)>bias_limit)
    			x_hat(2,0) = bias_limit;

    		if(x_hat(2,0)<-bias_limit)
    			x_hat(2,0) = -bias_limit;

    		alt.altitude = x_hat(0,0);
    		alt.vel_z = x_hat(1,0);

    		//test[0] = x_hat(0,0);

    		acc_z_pre = acc_z;
    		orb_publish(ORB_ID(alt_estimate), alt_estimate_pub, &alt);

    	}
    	// usleep(200000);		// 100ms
    	usleep(5000);
    }
    //关闭线程，串口
    warnx("[Alt_est] exiting.\n");

    thread_running = false;
    fflush(stdout);
    return 0;
}

#include <px4_config.h>
#include <px4_tasks.h>
#include <px4_posix.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <poll.h>
#include <string.h>

#include <uORB/uORB.h>
#include <uORB/topics/sensor_combined.h>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/myfirst_topic.h>

__EXPORT int myfirst_app_main(int argc,char *argv[]);

int myfirst_app_main(int argc,char *argv[])
{	
	orb_advert_t handle0;
	int handle1;
	struct myfirst_topic_s random0;
	struct myfirst_topic_s random1;
	handle0=orb_advertise(ORB_ID(myfirst_topic),&random0);
		//发布前将random中的数据进行公告
	PX4_INFO("uORB exercise");
	handle1=orb_subscribe(ORB_ID(myfirst_topic));//订阅自己发布的主题 (注意：一定要声明这个主题明，一般在对应的.h文件中声明)
	memset(&random0,0,sizeof(random0));
		random0.data1=rand()%1000;
		random0.data2=rand()%1000;
		random0.data3=rand()%1000;
	PX4_INFO("the new data is:\t%d\t%d\t%d",(int) random0.data1,(int) random0.data2,(int) random0.data3);
	//生成1000以内的随机数保存在结构体random0中
	orb_publish(ORB_ID(myfirst_topic),handle0,&random0);
		//发布主题
	
	PX4_INFO("the new data is:\t%d\t%d\t%d",(int) random0.data1,(int) random0.data2,(int) random0.data3);
	
	bool update;
	orb_check(handle1,&update);	
	
		//检查该主题的数据是否被获取过？
		//update=1,说明更新后未被获取
		//update=0,说明更新后已获取或者未更新
		if (update)
		{
			orb_copy(ORB_ID(myfirst_topic),handle1,&random1);
		 	PX4_INFO("the update data is:\t%d\t%d\t%d",(int) random0.data1,(int) random0.data2,(int) random0.data3);
		}
		else
		{
			PX4_WARN("date has not been updated");
		}
	return 0;
}

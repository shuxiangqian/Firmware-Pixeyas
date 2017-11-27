#include <px4_config.h>
#include <px4_tasks.h>
#include <px4_posix.h>
#include <unistd.h>
#include <stdio.h>

#include <stdlib.h>
#include <poll.h>
#include <string.h>

#include "stm32_gpio.h"
#include <termios.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <stdbool.h>
#include <errno.h>
#include <drivers/drv_hrt.h>
#include <uORB/topics/sonar_distance.h>

#define GPIO_UART8_TX_BREAK       (GPIO_OUTPUT|GPIO_OPENDRAIN|GPIO_SPEED_50MHz|GPIO_OUTPUT_CLEAR|GPIO_PORTE|GPIO_PIN1)


#define SONAR01                 0x01            // Address of the SFR01
#define SONAR02                 0x02            // Address of the SFR02
#define SONAR03                 0x03            // Address of the SFR03
#define SONAR04                 0x04            // Address of the SFR04
#define SONAR05                 0x05            // Address of the SFR05

#define RANGE_INCH 		0x50		//Real Ranging Mode - Result in inches
#define RANGE_CM		0x51		//Real Ranging Mode - Result in centimeters
#define QRANGE_INCH		0x53		//Real Ranging Mode - Result in inches, automatically Tx range back to controller as soon as ranging is complete.
#define QFASTRANGE_CM           0x54		//Real Ranging Mode - Result in centimeters, automatically Tx range back to controller as soon as ranging is complete.

#define FAKERANGE_INCH          0x56		//Fake Ranging Mode - Result in inches
#define FAKERANGE_CM            0x57		//Fake Ranging Mode - Result in centimeters
#define QFAKERANGE_INCH 	0x59		//Fake Ranging Mode - Result in inches, automatically Tx range back to controller as soon as ranging is complete.
#define QFAKERANGE_INCM         0x5A            //Fake Ranging Mode - Result in centimeters, automatically Tx range back to controller as soon as ranging is complete.

#define SENDBURST		0x5C		//Transmit an 8 cycle 40khz burst - no ranging takes place
#define GETVERSION 		0x5D            //Get software version - sends a single byte back to the controller
#define GETRANGE 		0x5E		//Get Range, returns two bytes (high byte first) from the most recent ranging.
#define GETSTATUS		0x5F		//Get Status, returns a single byte. Bit 0 indicates "Transducer locked", Bit 1 indicates "Advanced Mode"
#define SLEEP		 	0x60		//Sleep, shuts down the SRF01 so that it uses very low power (55uA).
#define UNLOCK			0x61		//Unlock. Causes the SRF01 to release and re-acquire its "Transducer lock". Used by our factory test routines.
#define ADVANCEDMODE            0x62		//Set Advanced Mode (Factory default) - When locked, SRF01 will range down to zero.
#define CLEARADVANCED           0x63		//Clear Advanced Mode - SRF01 will range down to approx. 12cm/5in,
#define BAUD19200		0x64		//Changes the baud rate to 19200
#define BAUD38400		0x65		//Changes the baud rate to 38400

#define FIRSTI2C		0xA0            //1st in sequence to change I2C address
#define THIRDI2C 		0xA5            //3rd in sequence to change I2C address
#define SECONDI2C 		0xAA            //2nd in sequence to change I2C address




//ORB_DEFINE(SonarGroupDistance, struct SonarGroupDistance_s,33,(const)__orb_SonarGroupDistance_fields);

static bool thread_should_exit = false;
static bool thread_running = false;
static int daemon_task;



// 声纳主函数
__EXPORT int sonargroup_service_main(int argc, char *argv[]);
int sonargroup_service_thread_main(int argc, char *argv[]);
//UART初始化，默认选择UART7
static int uart_init(char * uart_name);
// 设置波特率
static int set_uart_baudrate(const int fd, unsigned int baud);
// 发送SRF指令
static void SRF01(int fd_UART,unsigned char Address, unsigned char cmd);
// 读取测量距离
static int GetRange_new(int fd_UART,unsigned char Address);
// 使用提示函数
static void usage(const char *reason);

static void usage(const char *reason)
{
    if (reason) {
        fprintf(stderr, "%s\n", reason);
    }

    fprintf(stderr, "usage: sonargroup_service {start|stop|status} [param]\n\n");
    exit(1);
}


int sonargroup_service_main(int argc, char *argv[])
{
    if (argc < 2) {
        usage("[SONAR_GTOUP]missing command");
        return 1;
    }

    if (!strcmp(argv[1], "start")) {
        if (thread_running) {
            warnx("[SONAR_GTOUP]already running\n");
            return 0;
        }

        thread_should_exit = false;
        daemon_task = px4_task_spawn_cmd("sonar_group_service",
                                         SCHED_DEFAULT,
                                         SCHED_PRIORITY_MAX - 5,
                                         2000,
                                         sonargroup_service_thread_main,
                                         (argv) ? (char * const *)&argv[2] : (char * const *)NULL);
        return(0);
    }

    if (!strcmp(argv[1], "stop")) {
        thread_should_exit = true;
        return(0);
    }

    if (!strcmp(argv[1], "status")) {
        if (thread_running) {
            warnx("[SONAR_GTOUP]running");

        } else {
            warnx("[SONAR_GTOUP]stopped");
        }

        return(0);
    }

    usage("unrecognized command");
    return(1);
}

int sonargroup_service_thread_main(int argc, char *argv[])
{
    //初始化UARt8
    int uart_fd = uart_init("/dev/ttyS6");
    if(false == uart_fd)
    	return -1;
    if(false == set_uart_baudrate(uart_fd,9600)){
        warnx("[SONAR_GROUP]set_uart_baudrate is failed\n");
        return -1;
    }
    warnx("[SONAR_GROUP]UART init is successful\n");
    //static int count1=0;
    //线程启动标志
    thread_running = true;
    //初始化数据结构体
    struct sonar_distance_s sonar;
    memset(&sonar, 0, sizeof(sonar));
    sonar.count=0;
	sonar.distance_down[3]=0;
	sonar.distance_down[2]=0;
	sonar.distance_down[1]=0;
	sonar.distance_down[0]=0;
	sonar.vz=0.0f;
    //公告主题
    orb_advert_t SonarGroupDistance_pub = orb_advertise(ORB_ID(sonar_distance), &sonar);

    warnx("[SONAR_GROUP]service start successfully\n");
    SRF01(uart_fd,0,ADVANCEDMODE);
    usleep(1000);
    //获取所有超声波距离 Sonar_group GetRangeAll
    while(!thread_should_exit)
    {
        SRF01(uart_fd,2,RANGE_CM);
        //wait 70ms to recive the respons
        usleep(10000);
        SRF01(uart_fd,3,RANGE_CM);
		//wait 70ms to recive the respons
		usleep(10000);
		SRF01(uart_fd,4,RANGE_CM);
		//wait 70ms to recive the respons
		usleep(10000);
		SRF01(uart_fd,5,RANGE_CM);
		//wait 70ms to recive the respons
		usleep(10000);
		SRF01(uart_fd,1,RANGE_CM);
		//wait 70ms to recive the respons
		usleep(50000);

        sonar.count++;

		for (int i=1;i<=5;i++)
		{
			int range=GetRange_new(uart_fd,i);
			if (range!=-1979)
			{
				sonar.distance[i-1]=range;
				sonar.status[i-1]=1;
			}
			else
			{
				sonar.distance[i-1]=range;							//add by yly
				sonar.status[i-1]=0;
			}
		}

		sonar.distance_down[3] = sonar.distance_down[2];
		sonar.distance_down[2] = sonar.distance_down[1];
		sonar.distance_down[1] = sonar.distance_down[0];
		sonar.distance_down[0] = sonar.distance[0];

		if(sonar.count>=3)
		{
			if(sonar.distance_down[0]<0.0f)
			{
				sonar.distance_down[0] = sonar.distance_down[1];
				sonar.distance_filter = sonar.distance_down[0];
			}
			else
			{
				float k1=(sonar.distance_down[1]-sonar.distance_down[2])/100.0f;
				float k2=(sonar.distance_down[0]-sonar.distance_down[1])/100.0f;

				if(k2>0.2f&&k2>3*k1)
				{
					sonar.distance_filter = sonar.distance_down[1]+k1*0.2f;
					sonar.distance_down[0] = sonar.distance_filter;
				}
				else if(k2<-0.2f&&k2<3*k1)
				{
					sonar.distance_filter = sonar.distance_down[1]+k1*0.2f;
					sonar.distance_down[0] = sonar.distance_filter;
				}
				else
				{
					sonar.distance_filter = sonar.distance_down[0];
				}
			}

			sonar.vz=0.0f;
		}
		else
			sonar.distance_filter = sonar.distance[0];

		orb_publish(ORB_ID(sonar_distance), SonarGroupDistance_pub, &sonar);
    }
    //关闭线程，串口
    warnx("[SONAR_GROUP] exiting.\n");

    close(uart_fd);
    thread_running = false;
    fflush(stdout);
    return 0;
}


// Function to send commands to the SRF01
void SRF01(int fd_UART,unsigned char Address,unsigned char cmd)
{
  
  //配置成ＧＰＩＯ模式，发送高低电平“break"
    stm32_configgpio(GPIO_UART8_TX_BREAK);
    // Send a 2ms break to begin communications with the SRF01
    stm32_gpiowrite(GPIO_UART8_TX_BREAK,0);
    usleep(2000);
    stm32_gpiowrite(GPIO_UART8_TX_BREAK,1);
    usleep(10);

    //配置成URAT模式，发送指令
    stm32_configgpio(GPIO_UART8_TX);

  /*  ioctl(fd_UART, TIOCSBRK,0);
    usleep(2000000);
    ioctl(fd_UART, TIOCCBRK,0);
    //ioctl(fd_UART,TCSBRK,0);
*/
    write(fd_UART,&Address,sizeof(Address));
    write(fd_UART,&cmd,sizeof(cmd));

    // As RX and TX are the same pin it will have recieved the data we just sent out, as we dont want this we read it back and ignore it as junk before waiting for useful data to arrive
    usleep(3000);
    unsigned char buff2[20];
    read(fd_UART,&buff2,20);
}

int GetRange_new(int fd_UART,unsigned char Address)
{
    //发送测距指令
    SRF01(fd_UART,Address,GETRANGE);
    //等待1ms反馈
    //usleep(1000);//yly
    usleep(2000);

    //高低位两个字节
    unsigned char hByte,lByte;
    int readH=read(fd_UART,&hByte,1);
    int readL=read(fd_UART,&lByte,1);


    //如果接收不正常，反馈-1979错误代码
    if ((readH!=1) || (readL!=1))
    {
        return -1979;
    }else
    {
        int range = ((hByte<<8)+lByte);
        return range;
    }
}

int set_uart_baudrate(const int fd, unsigned int baud)
{
    int speed;

    switch (baud) {
    case 9600:   speed = B9600;   break;
    case 19200:  speed = B19200;  break;
    case 38400:  speed = B38400;  break;
    case 57600:  speed = B57600;  break;
    case 115200: speed = B115200; break;
    default:
        warnx("ERR: baudrate: %d\n", baud);
        return -EINVAL;
    }

    struct termios uart_config;

    int termios_state;

    /* fill the struct for the new configuration */
    tcgetattr(fd, &uart_config);

    /* clear ONLCR flag (which appends a CR for every LF) */
    uart_config.c_oflag &= ~ONLCR;
    /* no parity, one stop bit */
    uart_config.c_cflag &= ~(CSTOPB | PARENB);
    /* set baud rate */

    if ((termios_state = cfsetispeed(&uart_config, speed)) < 0) {
        warnx("ERR: %d (cfsetispeed)\n", termios_state);
        return false;
    }

    if ((termios_state = cfsetospeed(&uart_config, speed)) < 0) {
        warnx("ERR: %d (cfsetospeed)\n", termios_state);
        return false;
    }

    if ((termios_state = tcsetattr(fd, TCSANOW, &uart_config)) < 0) {
        warnx("ERR: %d (tcsetattr)\n", termios_state);
        return false;
    }
    /* Activate single wire mode */
    ioctl(fd, TIOCSSINGLEWIRE, SER_SINGLEWIRE_ENABLED);

    return true;
}

int uart_init(char * uart_name)
{
    // 可读写|如果欲打开的文件为终端机设备时，则不会将该终端机当成进程控制终端机|非阻塞
    int serial_fd = open(uart_name, O_RDWR | O_NOCTTY | O_NONBLOCK);

    if (serial_fd < 0) {
        err(1, "failed to open port: %s", uart_name);
        return false;
    }
    return serial_fd;
}

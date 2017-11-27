/****************************************************************************
 *
 *   Copyright (c) 2012-2016 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file px4_simple_app.c
 * Minimal application example for PX4 autopilot
 *
 * @author Example User <mail@example.com>
 */

#include <px4_config.h>
#include <px4_tasks.h>
#include <px4_posix.h>
#include <unistd.h>
#include <stdio.h>
#include <poll.h>
#include <string.h>

#include <uORB/uORB.h>
#include <uORB/topics/sensor_combined.h>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/sonar_distance.h>
#include <uORB/topics/manual_control_setpoint.h>
//#include <uORB/topics/alt_ctrl.h>

__EXPORT int px4_simple_app_main(int argc, char *argv[]);

int px4_simple_app_main(int argc, char *argv[])
{
	PX4_INFO("Hello Sky!");

	/* subscribe to sensor_combined topic */
	//int sensor_sub_fd = orb_subscribe(ORB_ID(sensor_combined));

	/* advertise attitude topic */
	//struct vehicle_attitude_s att;
	//memset(&att, 0, sizeof(att));
	//orb_advert_t att_pub = orb_advertise(ORB_ID(vehicle_attitude), &att);
	int sensor_sub_fd = orb_subscribe(ORB_ID(sonar_distance));
	orb_set_interval(sensor_sub_fd, 50);	// limit the update rate to 5 Hz
	struct sonar_distance_s sonar;
	memset(&sonar, 0, sizeof(sonar));

	//int alt_ctrl_sub_ = orb_subscribe(ORB_ID(alt_ctrl));
	//struct alt_ctrl_s alt_control;
	//memset(&alt_control,0,sizeof(alt_control));

	/*int sensor_sub_fd=orb_subscribe(ORB_ID(manual_control_setpoint));
	struct manual_control_setpoint_s manual;
	memset(&manual, 0, sizeof(manual));	*/
	/* one could wait for multiple topics with this technique, just using one here */
	px4_pollfd_struct_t fds[] = {
		{ .fd = sensor_sub_fd,   .events = POLLIN },
		/* there could be more file descriptors here, in the form like:
		 * { .fd = other_sub_fd,   .events = POLLIN },
		 */
	};

	int error_counter = 0;

	for (int i = 0; i < 1000; i++)
	{
		/* wait for sensor update of 1 file descriptor for 1000 ms (1 second) */
		int poll_ret = px4_poll(fds, 1, 70);
		/* handle the poll result */
		if (poll_ret == 0)
		{
			/* this means none of our providers is giving us data */
			PX4_ERR("Got no data within a second");

		}
		else if (poll_ret < 0)
		{
			/* this is seriously bad - should be an emergency */
			if (error_counter < 10 || error_counter % 50 == 0) {
				/* use a counter to prevent flooding (and slowing us down) */
				PX4_ERR("ERROR return value from poll(): %d", poll_ret);
			}
			error_counter++;
		}
		else
		{
			if (fds[0].revents & POLLIN)
			{
				/* obtained data for the first file descriptor */
				//struct sensor_combined_s raw;
				//struct sonar_distance_s raw;
				/* copy sensors raw data into local buffer */
				orb_copy(ORB_ID(sonar_distance), sensor_sub_fd, &sonar);
				//orb_copy(ORB_ID(alt_ctrl), alt_ctrl_sub_, &alt_control);

				//printf("----------------------------------\n");
				for(int j=1;j<=5;j++)
				{
					printf("[SONAR]Range(%d)=%d(cm)\n",j,sonar.distance[j-1]);
					printf("[SONAR]status(%d)=%d\n\n",j,sonar.status[j-1]);
					//warnx("[SD] thrust=%.2f alt_sp=%.2f alt_now=%.2f",(double)alt_control.thrust,(double)alt_control.alt_sp,(double)alt_control.alt_measure);

				}
				printf("=========Press CTRL+C to abort=========\n");

				//orb_copy(ORB_ID(manual_control_setpoint),sensor_sub_fd,&manual);
				//warnx("manual:x=%.2f,y=%.2f,z=%.2f",(double)manual.x,(double)manual.y,(double)manual.z);
				//warnx("=========Press CTRL+C to abort=========\n");

				//orb_publish(ORB_ID(vehicle_attitude), att_pub, &att);
			}
			/* there could be more file descriptors here, in the form like:
			 * if (fds[1..n].revents & POLLIN) {}
			 */
		}
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

	PX4_INFO("exiting");

	return 0;
}

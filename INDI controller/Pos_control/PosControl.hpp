/****************************************************************************
 *
 *   Copyright (c) 2018 - 2019 PX4 Development Team. All rights reserved.
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
 * @file PosControl.hpp
 *Author:Huazi Cao
 */

#pragma once
#include "ControlMath.hpp"
#include <string.h>
#include <lib/mathlib/mathlib.h>
#include <matrix/matrix/math.hpp>
#include <uORB/uORB.h>
#include <uORB/Publication.hpp>
#include <uORB/Subscription.hpp>
#include <uORB/topics/motor_speed.h>
#include <uORB/topics/vehicle_attitude_setpoint.h>
#include <uORB/topics/vehicle_local_position_setpoint.h>
#include <uORB/topics/pos_helper.h> // load log
#include <float.h>
#include <px4_platform_common/defines.h>
#include <geo/geo.h>


#define G	        9.8066f
#define SUM_max_x       0.85f                                          //x方向水平最大积分
#define SUM_max_y       0.85f
#define SUM_max_z       0.85f
#define Pi              3.1415926f
#define or_rotorcraft_max_rotors 4

struct PositionControlStates {
        matrix::Vector3f position;
        matrix::Vector3f velocity;
        matrix::Vector3f acceleration;
        float yaw;
};

// User define control parameters
struct ControlParas{                                               //控制器
     matrix::Matrix<float,3, 3> K_x;
     matrix::Matrix<float,3, 3> K_v;
     matrix::Matrix<float,3, 3> K_a;
};

// User define structure.
struct Autopilot{
     matrix::Vector3f pos_err;                                    //位置误差
     matrix::Vector3f vel_err;                                    //速度误差
     matrix::Vector3f acc_err;                                    //加速度误差
     matrix::Vector3f a_c;                                        //加速度控制量
     matrix::Vector3f tau_bz_f;
     matrix::Vector3f tau_bz_c;                                   //拉力控制量
     matrix::Vector3f f_iusl;                                     //输出推力
};

// User define structure.
struct ema_filter_s{
    double scale[9];
    double in[3];
    double alpha[3];
    double out[3];
};

// User define structure.
struct ema_filter_generic_s{
    double scale[8*8];
    double in[8];
    double alpha[8];
    double out[8];
};

class PosControl
{
public:

        PosControl() = default;
        ~PosControl() = default;

        /**
         * Set the parameters about propeller to force and torque
         * get the parameters from test data
         */
        void setpropParams(const float prop2force, const float prop2torque){_K_f = prop2force; _K_t = prop2torque;};
        /**
         * Set the controller gains and ESO gains
         * Find detials from my paper
         */
        void setMasses(const float m_multrotor, const float m_manipulator){_Mb = m_multrotor; _SUM_mi = m_manipulator;};
        /**
         * Set the hover thrust of the uav
         * Find detials from my paper
         */
        void setUerDeineHoverThrust(const float thr_hover){_hover_thr = thr_hover;};

        void setControlParas(const matrix::Vector3f &Kx, const matrix::Vector3f &Kv, const matrix::Vector3f &Ka);

        void setESOParas(const matrix::Vector3f &ESO_v);
        /**
         * Set the maximum velocity to execute with feed forward and position control
         * @param vel_horizontal horizontal velocity limit
         * @param vel_up upwards velocity limit
         * @param vel_down downwards velocity limit
         */
        void setVelocityLimits(const float vel_horizontal, const float vel_up, float vel_down);

        /**
         * Set the minimum and maximum collective normalized thrust [0,1] that can be output by the controller
         * @param min minimum thrust e.g. 0.1 or 0
         * @param max maximum thrust e.g. 0.9 or 1
         */
        void setThrustLimits(const float min, const float max);

        /**
         * Set the maximum tilt angle in radians the output attitude is allowed to have
         * @param tilt angle in radians from level orientation
         */
        void setTiltLimit(const float tilt) { _lim_tilt = tilt; }

        /**
         * Set the normalized hover thrust, from hover estimate
         * @param thrust [0.1, 0.9] with which the vehicle hovers not acelerating down or up with level orientation
         */
        void setHoverThrust(const float hover_thrust) { _hover_thrust = math::constrain(hover_thrust, 0.1f, 0.9f); }

        /**
         * Update the hover thrust without immediately affecting the output
         * by adjusting the integrator. This prevents propagating the dynamics
         * of the hover thrust signal directly to the output of the controller.
         */
        void updateHoverThrust(const float hover_thrust_new);

        /**
         * Pass the current vehicle state to the controller
         * @param PositionControlStates structure
         */
        void setState(const PositionControlStates &states);

        /**
         * Pass the current static estimate
         */
        void setStaticDisturbance(const matrix::Vector3f &static_af);

        /**
         * Pass the desired setpoints
         * Note: NAN value means no feed forward/leave state uncontrolled if there's no higher order setpoint.
         * @param setpoint a vehicle_local_position_setpoint_s structure
         */
        void setInputSetpoint(const vehicle_local_position_setpoint_s &setpoint);

        /**
         * set EMA filter parameters
         * Note: This function is only used for 3-dimensional vectors
         */
        void setEMAFilterParams();

        /**
         * set EMA generic filter parameters
         * Note: This function is used for any-dimensional vectors
         */
        void setEMAFilterGenericParams();

        /**
         * get motor speed
         */
        void setMotorSpeed(const motor_speed_s &motor_speed);

        /**
         * Apply P-position and PID-velocity controller that updates the member
         * thrust, yaw- and yawspeed-setpoints.
         * @see _thr_sp
         * @see _yaw_sp
         * @see _yawspeed_sp
         * @param dt time in seconds since last iteration
         * @return true if update succeeded and output setpoint is executable, false if not
         */
        bool update(const matrix::Quatf &q, const matrix::Vector3f &acc, const float dt);

        /**
         * Get the controllers output local position setpoint
         * These setpoints are the ones which were executed on including PID output and feed-forward.
         * The acceleration or thrust setpoints can be used for attitude control.
         * @param local_position_setpoint reference to struct to fill up
         */
        void getLocalPositionSetpoint(vehicle_local_position_setpoint_s &local_position_setpoint) const;

        /**
         * Get the controllers output attitude setpoint
         * This attitude setpoint was generated from the resulting acceleration setpoint after position and velocity control.
         * It needs to be executed by the attitude controller to achieve velocity and position tracking.
         * @param attitude_setpoint reference to struct to fill up
         */
        void getAttitudeSetpoint(vehicle_attitude_setpoint_s &attitude_setpoint) const;

        /**
         * Exponential moving average filter(EMA)
         * @param raw Raw data that needs to be filtered
         * @param scale
         * @param in
         * @param alpha
         * @param out data after filtering
         */
        void ema_filter(const double raw[3], const double scale[3*3], double in[3], const double alpha[3], double out[3]);

        /**
         * Exponential moving average filter(EMA)
         * @param raw Raw data that needs to be filtered
         * @param raw_size size of the raw data
         * @param scale
         * @param scale_rows
         * @param scale_cols
         * @param in
         * @param alpha
         * @param out data after filtering
         */
        void ema_filter_generic(const double raw[], size_t raw_size,
                               const double scale[], size_t scale_rows, size_t scale_cols,
                               double in[], const double alpha[], double out[]);
private:
        bool _updateSuccessful();

        void _positionControl(const matrix::Quatf &q, const matrix::Vector3f &acc, const float dt); ///< Position control. Details can be found in my paper.

        // Gains
        ControlParas _contolParas;

        // lpf parameters
        ema_filter_s _ema_filter;
        ema_filter_generic_s _ema_generic_filter;

        //uORB data for logger
        uORB::Publication<pos_helper_s>	pos_helper_pub{ORB_ID(pos_helper)};     //输出uorb

        // User defined varibles
        Autopilot _autopilot;

        // parameters about propellers
        float _K_f;
        float _K_t;
        // Masses of multirotor and manipulator
        float _Mb;
        float _SUM_mi;
        float _hover_thr;
        matrix::Vector3f g = matrix::Vector3f(0, 0,G);

        // User defined log
        pos_helper_s pos_helper {};

        motor_speed_s _motor_speed;
        // Limits
        float _lim_vel_horizontal{}; ///< Horizontal velocity limit with feed forward and position control
        float _lim_vel_up{}; ///< Upwards velocity limit with feed forward and position control
        float _lim_vel_down{}; ///< Downwards velocity limit with feed forward and position control
        float _lim_thr_min{}; ///< Minimum collective thrust allowed as output [-1,0] e.g. -0.9
        float _lim_thr_max{}; ///< Maximum collective thrust allowed as output [-1,0] e.g. -0.1
        float _lim_tilt{}; ///< Maximum tilt from level the output attitude is allowed to have

        float _hover_thrust{}; ///< Thrust [0.1, 0.9] with which the vehicle hovers not accelerating down or up with level orientation

        // States
        matrix::Vector3f _pos; /**< current position */
        matrix::Vector3f _vel; /**< current velocity */
        matrix::Vector3f _acc; /**< current acceleration */
        matrix::Vector3f _vel_dot; /**< velocity derivative (replacement for acceleration estimate) */
        matrix::Vector3f _vel_int; /**< integral term of the velocity controller */
        float _yaw{}; /**< current heading */
        matrix::Vector3f _static_af = matrix::Vector3f(0, 0,0);

        // Setpoints
        matrix::Vector3f _pos_sp; /**< desired position */
        matrix::Vector3f _vel_sp; /**< desired velocity */
        matrix::Vector3f _acc_sp; /**< desired acceleration */
        matrix::Vector3f _acc_cal; /**< desired acceleration */
        matrix::Vector3f _thr_sp; /**< desired thrust */
        float _yaw_sp{}; /**< desired heading */
        float _yawspeed_sp{}; /** desired yaw-speed */
        matrix::Vector3f delta_v;
        struct usr_ESO
        {
        matrix::Vector3f pos_est;
        matrix::Vector3f vel_est;
        matrix::Vector3f delta_est;
        matrix::Matrix3f gain_ESO;
        void setZero(matrix::Vector3f pos,matrix::Vector3f vel)
        {
                pos_est = pos;
                vel_est  = vel;
                delta_est.setZero();
        }
        } _usr_eso;
        void PositionESO(matrix::Vector3f pos_in,matrix::Vector3f f,float dt);
};

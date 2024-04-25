#pragma once

#include <algorithm>
#include <math.h>
#include <mathlib/mathlib.h>
#include <matrix/matrix/math.hpp>
#include <mathlib/math/Limits.hpp>
#include <lib/mixer/MultirotorMixer/MultirotorMixer.hpp>
#include <uORB/topics/rate_ctrl_status.h>
#include <uORB/topics/motor_speed.h>


//log
#include <uORB/uORB.h>
#include <uORB/Publication.hpp>
#include <uORB/topics/att_helper.h> // load log

#define or_rotorcraft_max_rotors 4


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

class Att_Control
{
public:
        Att_Control() = default;
        ~Att_Control() = default;

        void Controllerinit();

        /**
         * Set attitude controller and attitude ESO gains
         */
        void setControllerGain(const matrix::Vector3f &KQ,const matrix::Vector3f &KW);
        /**
         * Set attitude controller and attitude ESO gains
         */
        void setInertiaMatrix(const matrix::SquareMatrix<float, 3> &Ib){_I_b = Ib;_I_b_inve = matrix::inv(Ib);}

        /**
         * Set the parameters about propeller to force and torque
         * get the parameters from test data
         */
        void setpropParams(const float prop2force, const float prop2torque, const double arm_length);

        /**
         * Set hard limit for output rate setpoints
         * @param rate_limit [rad/s] 3D vector containing limits for roll, pitch, yaw
         */
        void setRateLimit(const matrix::Vector3f &rate_limit) { _rate_limit = rate_limit; }
        /**
         * Set the integral item of the ESO as zero
         */
        void resetESO();

        /**
         * Set a new attitude setpoint replacing the one tracked before
         * @param qd desired vehicle attitude setpoint
         * @param yawspeed_setpoint [rad/s] yaw feed forward angular rate in world frame
         */
        void setAttitudeSetpoint(const matrix::Quatf &qd, const float yawspeed_setpoint) { _attitude_setpoint_q = qd; _attitude_setpoint_q.normalize(); _yawspeed_setpoint = yawspeed_setpoint; }

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
         * Adjust last known attitude setpoint by a delta rotation
         * Optional use to avoid glitches when attitude estimate reference e.g. heading changes.
         * @param q_delta delta rotation to apply
         */
        void adaptAttitudeSetpoint(const matrix::Quatf &q_delta) { _attitude_setpoint_q = q_delta * _attitude_setpoint_q; }

        /**
         * Run one control loop cycle calculation
         * @param q estimation of the current vehicle attitude unit quaternion
         * @return [rad/s] body frame 3D angular rate setpoint vector to be executed by the rate controller
         */
        void update(const matrix::Quatf &q, const matrix::Vector3f &rate, const matrix::Vector3f &tau_static,
                          const float &dt, const bool &landed, matrix::Vector3f &torque,matrix::Vector3f &rates_sp);
        /**
         * Set saturation status
         * @param status message from mixer reporting about saturation
         */
        void getRateControlStatus(rate_ctrl_status_s &rate_ctrl_status);

        /**
         * get motor speed
         */
        void setMotorSpeed(const motor_speed_s &motor_speed);

        /**
         * Exponential moving average filter(EMA)
         * @param raw Raw data that needs to be filtered
         * @param scale
         * @param in
         * @param alpha
         * @param out data after filtering
         */
        static void ema_filter(const double raw[3], const double scale[3*3], double in[3], const double alpha[3], double out[3]);

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
        static void ema_filter_generic(const double raw[], size_t raw_size,
                               const double scale[], size_t scale_rows, size_t scale_cols,
                               double in[], const double alpha[], double out[]);
private:

        /**
         * Set hard limit for output rate setpoints
         * @param rate_limit [rad/s] 3D vector containing limits for roll, pitch, yaw
         */
        void runAttitudeControl(const matrix::Quatf &q, const matrix::Vector3f &rate,const matrix::Vector3f &tau_static,
                                const float &dt, matrix::Vector3f& torque,matrix::Vector3f &rates_sp);
        /**
         * Run attitude ESO
         * @param angular rate[rad/s] output tau[N*M]
         */
        void UsrAttitudeESO(matrix::Vector3f bm_omega,matrix::Vector3f u,float dt);

        matrix::Vector3f _rate_limit;
        float _yaw_w{0.f}; ///< yaw weight [0,1] to deprioritize caompared to roll and pitch

        matrix::Quatf _attitude_setpoint_q; ///< latest known attitude setpoint e.g. from position control
        float _yawspeed_setpoint{0.f}; ///< latest known yawspeed feed-forward setpoint

        matrix::Matrix<float,3, or_rotorcraft_max_rotors> _R_Torque;
        //uORB data for logger
        uORB::Publication<att_helper_s>	att_helper_pub{ORB_ID(att_helper)};     //输出uorb

        // lpf parameters
        ema_filter_s _ema_filter;
        ema_filter_generic_s _ema_generic_filter;
        motor_speed_s _motor_speed;

        matrix::Vector3f _wf;
        matrix::Vector3f _wf_prev;
        matrix::Vector3f _dot_Omega_f;

        struct usr_ESO
        {
        matrix::Vector3f bm_rate_esti;                             //ESO
        matrix::Vector3f bm_rate_esti_dot;
        matrix::Vector3f bm_gain;
        matrix::Vector3f delta_esti;
        matrix::Vector3f delta_esti_dot;
        } _usr_eso;


        struct usr_att_controller                        //controller
        {
        matrix::Matrix<float,3, 3> K_xi;
        matrix::Matrix<float,3, 3> K_omega;
        } _controller_param;

        matrix::Vector3f _tau;
        matrix::Matrix<float,3, 3> _I_b;
        matrix::SquareMatrix<float, 3> _I_b_inve;


        // User defined log
        att_helper_s _att_helper {};

};

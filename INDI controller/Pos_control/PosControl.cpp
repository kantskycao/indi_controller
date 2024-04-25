/*****************************************************************************/

/**
* @file PosControl.cpp
* Author:Huazi Cao
*/

#include "PosControl.hpp"
#include <mathlib/mathlib.h>

using namespace matrix;
using namespace time_literals;

void PosControl::setControlParas(const matrix::Vector3f &Kx, const matrix::Vector3f &Kv, const matrix::Vector3f &Ka)
{
    _contolParas.K_x(0,0) = Kx(0);
    _contolParas.K_x(1,1) = Kx(1);
    _contolParas.K_x(2,2) = Kx(2);
    _contolParas.K_v(0,0) = Kv(0);
    _contolParas.K_v(1,1) = Kv(1);
    _contolParas.K_v(2,2) = Kv(2);
    _contolParas.K_a(0,0) = Ka(0);
    _contolParas.K_a(1,1) = Ka(1);
    _contolParas.K_a(2,2) = Ka(2);
}
void PosControl::setESOParas(const matrix::Vector3f &ESO_v)
{
    _usr_eso.gain_ESO(0,0)= ESO_v(0);
    _usr_eso.gain_ESO(1,1)= ESO_v(1);
    _usr_eso.gain_ESO(2,2)= ESO_v(2);
}

void PosControl::setVelocityLimits(const float vel_horizontal, const float vel_up, const float vel_down)
{
        _lim_vel_horizontal = vel_horizontal;
        _lim_vel_up = vel_up;
        _lim_vel_down = vel_down;
}

void PosControl::setThrustLimits(const float min, const float max)
{
        // make sure there's always enough thrust vector length to infer the attitude
        _lim_thr_min = math::max(min, 10e-4f);
        _lim_thr_max = max;
}

void PosControl::updateHoverThrust(const float hover_thrust_new)
{
        _vel_int(2) += (hover_thrust_new - _hover_thrust) * (CONSTANTS_ONE_G / hover_thrust_new);
        setHoverThrust(hover_thrust_new);
}

void PosControl::setState(const PositionControlStates &states)
{
        _pos = states.position;
        _vel = states.velocity;
        _acc = states.acceleration;
        _yaw = states.yaw;
        _vel_dot = states.acceleration;
}

void PosControl::setStaticDisturbance(const matrix::Vector3f &static_af)
{
        _static_af = static_af;
}

void PosControl::setInputSetpoint(const vehicle_local_position_setpoint_s &setpoint)
{
        _pos_sp = Vector3f(setpoint.x, setpoint.y, setpoint.z);
        _vel_sp = Vector3f(setpoint.vx, setpoint.vy, setpoint.vz);
        _acc_sp = Vector3f(setpoint.acceleration);
        _yaw_sp = setpoint.yaw;
        _yawspeed_sp = setpoint.yawspeed;
}
void PosControl::setMotorSpeed(const motor_speed_s &motor_speed)
{
    _motor_speed = motor_speed;
}

void PosControl::setEMAFilterParams(){
    double scale_temp[9] = {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0};
    double in_temp[3] = { 0.0, 0.0, 0.0 };
    double alpha_temp[3] = { 0.001, 0.001, 0.001 };
    double out_temp[3] = { 0.0, 0.0, 0.0 };

    std::copy(std::begin(scale_temp), std::end(scale_temp), _ema_filter.scale);
    std::copy(std::begin(in_temp), std::end(in_temp), _ema_filter.in);
    std::copy(std::begin(alpha_temp), std::end(alpha_temp), _ema_filter.alpha);
    std::copy(std::begin(out_temp), std::end(out_temp), _ema_filter.out);

}
void PosControl::setEMAFilterGenericParams(){
    double scale_temp[8*8] = {
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    double in_temp[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0 };
    double alpha_temp[8] = { 0.001, 0.001, 0.001 };
    double out_temp[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0 };

    std::copy(std::begin(scale_temp), std::end(scale_temp), _ema_generic_filter.scale);
    std::copy(std::begin(in_temp), std::end(in_temp), _ema_generic_filter.in);
    std::copy(std::begin(alpha_temp), std::end(alpha_temp), _ema_generic_filter.alpha);
    std::copy(std::begin(out_temp), std::end(out_temp), _ema_generic_filter.out);
}

/* --- Filter Definition ------------------------------------------------------- */
/* Exponential moving average filter(EMA) generic */
void PosControl::ema_filter(const double raw[3], const double scale[3*3],
           double in[3], const double alpha[3], double out[3])
{
  double v[3];
  int i;

  for(i = 0; i < 3; i++)
    v[i] = raw[i];

  for(i = 0; i < 3; i++) {
    in[i] =
      scale[3*i + 0] * v[0] +
      scale[3*i + 1] * v[1] +
      scale[3*i + 2] * v[2];

    if (!std::isnan(out[i]))
      out[i] += alpha[i] * (in[i] - out[i]);
    else
      out[i] = in[i];
  }
}

void PosControl::ema_filter_generic(const double raw[], size_t raw_size,
                        const double scale[], size_t scale_rows, size_t scale_cols,
                        double in[], const double alpha[], double out[])
{
  // size_t raw_size = sizeof(raw) / sizeof(raw[0]);
  // size_t scale_rows = sizeof(scale) / (sizeof(scale[0]) * raw_size);
  // size_t scale_cols = raw_size;
  double v[raw_size];
  size_t i, j;

  for (i = 0; i < raw_size; i++)
    v[i] = raw[i];

  for (i = 0; i < raw_size; i++) {
    in[i] = 0.0;
    for (j = 0; j < scale_cols; j++) {
      in[i] += scale[scale_cols * i + j] * v[j];
    }

    if (!std::isnan(out[i]))
      out[i] += alpha[i] * (in[i] - out[i]);
    else
      out[i] = in[i];
  }
}

// TODO:I will rewrite this function
bool PosControl::update(const matrix::Quatf &q, const matrix::Vector3f &acc, const float dt)
{
        // x and y input setpoints always have to come in pairs
        const bool valid = (PX4_ISFINITE(_pos_sp(0)) == PX4_ISFINITE(_pos_sp(1)))
                           && (PX4_ISFINITE(_vel_sp(0)) == PX4_ISFINITE(_vel_sp(1)))
                           && (PX4_ISFINITE(_acc_sp(0)) == PX4_ISFINITE(_acc_sp(1)));

        // if(!PX4_ISFINITE(_pos_sp(0))&&!PX4_ISFINITE(_pos_sp(1))){
	// 	if(PX4_ISFINITE(_vel_sp(0)) && PX4_ISFINITE(_vel_sp(1))){
	// 		if (PX4_ISFINITE(_pos(0))&&PX4_ISFINITE(_pos(2))){
	// 			_pos_sp =  _pos + _vel_sp*dt;
	// 		}
	// 	}
	// }
        // if only velocity cmd is inputed, we integret velocity cmd to obtain position cmd.
        // This procedure is necessary, because the whole position controller requires position cmd.
        // unless noly velocity control is required by user
       for(int i = 0; i < 3; i++){
                if ((!PX4_ISFINITE(_pos_sp(i))) && PX4_ISFINITE(_vel_sp(i)) &&PX4_ISFINITE(_pos(i))) {
                        // pos_sp is NAN, vel_sp is not NAN
                        _pos_sp(i) = _pos(i) + _vel_sp(i)*dt;
                }
       }
        _positionControl(q, acc, dt);

        _yawspeed_sp = PX4_ISFINITE(_yawspeed_sp) ? _yawspeed_sp : 0.f;
        _yaw_sp = PX4_ISFINITE(_yaw_sp) ? _yaw_sp : _yaw; // TODO: better way to disable yaw control

//        PX4_WARN("valid = %d\n", valid);
//        PX4_WARN("_updateSuccessful = %d\n", _updateSuccessful());
        return valid && _updateSuccessful();
}

void PosControl::_positionControl(const matrix::Quatf &q,const matrix::Vector3f &acc, const float dt)
{
        // calculating R matrix from body to world frame
        Dcmf Rb2w(q);

        matrix::Vector3f a_f = Rb2w * acc + Vector3f(0, 0, G);
        double a_arr[3] = {a_f(0), a_f(1), a_f(2)};

        ema_filter(a_arr, _ema_filter.scale, _ema_filter.in, _ema_filter.alpha, _ema_filter.out);
        Vector3f acc_f{ float(_ema_filter.out[0]), float(_ema_filter.out[1]), float(_ema_filter.out[2])};

        // calculating tracking error
        _autopilot.pos_err       = _pos_sp - _pos;
        // make sure there are no NAN elements for further reference while constraining
        ControlMath::setZeroIfNanVector3f(_autopilot.pos_err);
        // check _vel_sp
        ControlMath::setZeroIfNanVector3f(_vel_sp);
        // check _acc_sp
        ControlMath::setZeroIfNanVector3f(_acc_sp);
        // calculating tracking error
        _autopilot.vel_err       = _vel_sp - _vel;
        _autopilot.acc_err       = _acc_sp - acc_f;
        //此处的加速度值来源于local_position的微分项，大概率不准确，应该使用IMU加速度计。

        //外扰数据_static_af
        ControlMath::setZeroIfNanVector3f(_static_af);
        _static_af(0) = math::constrain(_static_af(0), -5.0f, 5.0f);
        _static_af(1) = math::constrain(_static_af(1), -5.0f, 5.0f);
        _static_af(2) = math::constrain(_static_af(2), -25.0f, 25.0f);

        _autopilot.a_c = _contolParas.K_x * _autopilot.pos_err + _contolParas.K_v * _autopilot.vel_err + _contolParas.K_a * _autopilot.acc_err + _acc_sp;

        double wprop2_arr[or_rotorcraft_max_rotors];
        for (int i = 0; i < or_rotorcraft_max_rotors; ++i) {
            wprop2_arr[i] = _motor_speed.motor_spd_rpm[i] * _motor_speed.motor_spd_rpm[i];
        }
        ema_filter_generic(wprop2_arr, or_rotorcraft_max_rotors, _ema_generic_filter.scale,  or_rotorcraft_max_rotors,
                            or_rotorcraft_max_rotors, _ema_generic_filter.in, _ema_generic_filter.alpha, _ema_generic_filter.out);        //电机转速进行滤波

        double T_bz_f = 0.0;
        for (int i = 0; i < or_rotorcraft_max_rotors; ++i)
            T_bz_f += _K_f * _ema_generic_filter.out[i] * _ema_generic_filter.out[i];

        _autopilot.tau_bz_f = Rb2w *  Vector3f(0.0, 0.0, T_bz_f / _Mb);
        _autopilot.tau_bz_c = _autopilot.tau_bz_f + _autopilot.a_c - acc_f;

        _autopilot.f_iusl = -_Mb *  _autopilot.tau_bz_c;

        //convert force to acceleration
        _acc_cal = -_autopilot.f_iusl/(_Mb) + g;



        // Assume standard acceleration due to gravity in vertical direction for attitude generation
        Vector3f body_z = Vector3f(-_acc_cal(0), -_acc_cal(1), CONSTANTS_ONE_G).normalized();

        ControlMath::limitTilt(body_z, Vector3f(0, 0, 1), _lim_tilt);

        /* Scale thrust assuming hover thrust produces standard gravity
        * from PX4:        float collective_thrust = _acc_sp(2) * (_hover_thrust / CONSTANTS_ONE_G) - _hover_thrust;
        * old version:     float collective_thrust = _acc_sp(2) * (_hover_thrust / G) - _hover_thrust;
        * newest version:  float collective_thrust = (_acc_sp(2) - G) * (_Mb + _SUM_mi);
        * result: PX4 is right. Only the codes from PX4 can run well.
        */
        float collective_thrust = _acc_cal(2) * (_hover_thr / CONSTANTS_ONE_G) - _hover_thr;
        // float collective_thrust = _acc_cal(2) *(_hover_thrust / CONSTANTS_ONE_G) - _hover_thrust;
        // Project thrust to planned body attitude
        collective_thrust /= (Vector3f(0, 0, 1).dot(body_z));
        collective_thrust = math::min(collective_thrust, -_lim_thr_min);
        _thr_sp = body_z * collective_thrust;

        // Saturate maximal vertical thrust
        _thr_sp(2) = math::max(_thr_sp(2) , -_lim_thr_max);

        // Get allowed horizontal thrust after prioritizing vertical control
        const float thrust_max_squared = _lim_thr_max * _lim_thr_max;
        const float thrust_z_squared = _thr_sp(2) * _thr_sp(2) ;
        const float thrust_max_xy_squared = thrust_max_squared - thrust_z_squared;
        float thrust_max_xy = 0;

        if (thrust_max_xy_squared > 0) {
                thrust_max_xy = sqrtf(thrust_max_xy_squared);
        }

        // Saturate thrust in horizontal direction
        matrix::Vector2f thrust_sp_xy(_thr_sp(0),_thr_sp(1));

        float thrust_sp_xy_norm = thrust_sp_xy.norm();

        //Vector3f _thr_sp;
        if (thrust_sp_xy_norm > thrust_max_xy) {
                _thr_sp(0) = thrust_sp_xy(0) / thrust_sp_xy_norm * thrust_max_xy;
                _thr_sp(1) = thrust_sp_xy(1) / thrust_sp_xy_norm * thrust_max_xy;
        }
        // Integrator anti-windup in vertical direction
	// if ((_thr_sp(2) >= -_lim_thr_min && _autopilot.pos_err(2) <= 0.0f) ||
	//     (_thr_sp(2) <= -_lim_thr_max && _autopilot.pos_err(2) >= 0.0f)) {
	// 	_autopilot.pos_err(2) = 0.f;
	// }

        // pos_helper_pub.publish(pos_helper);
}

bool PosControl::_updateSuccessful()
{
        bool valid = true;

        // For each controlled state the estimate has to be valid
        for (int i = 0; i <= 2; i++) {
                if (PX4_ISFINITE(_pos_sp(i))) {
                        valid = valid && PX4_ISFINITE(_pos(i));
                }

                if (PX4_ISFINITE(_vel_sp(i))) {
                        valid = valid && PX4_ISFINITE(_vel(i)) && PX4_ISFINITE(_vel_dot(i));
                }
        }

        // There has to be a valid output accleration and thrust setpoint otherwise there was no
        // setpoint-state pair for each axis that can get controlled
        valid = valid && PX4_ISFINITE(_acc_cal(0)) && PX4_ISFINITE(_acc_cal(1)) && PX4_ISFINITE(_acc_cal(2));
        valid = valid && PX4_ISFINITE(_thr_sp(0)) && PX4_ISFINITE(_thr_sp(1)) && PX4_ISFINITE(_thr_sp(2));
        return valid;
}

void PosControl::PositionESO(matrix::Vector3f pos_in,matrix::Vector3f u,float dt)
{
    // check
    ControlMath::setZeroIfNanVector3f(_usr_eso.pos_est);
    ControlMath::setZeroIfNanVector3f(_usr_eso.vel_est);
    ControlMath::setZeroIfNanVector3f(_usr_eso.delta_est);
    ControlMath::setZeroIfNanVector3f(pos_in);
    ControlMath::setZeroIfNanVector3f(u);
    dt = PX4_ISFINITE(dt)? dt:0.0f;
    Vector3f p_est_dot = _usr_eso.vel_est + 3.0f * _usr_eso.gain_ESO * (pos_in - _usr_eso.pos_est);
    Vector3f v_est_dot = u + _usr_eso.delta_est  + 3.0f * _usr_eso.gain_ESO * _usr_eso.gain_ESO * (pos_in - _usr_eso.pos_est);
    Vector3f d_est_dot = _usr_eso.gain_ESO * _usr_eso.gain_ESO * _usr_eso.gain_ESO * (pos_in - _usr_eso.pos_est);
    _usr_eso.pos_est = _usr_eso.pos_est  + p_est_dot*dt;
    _usr_eso.vel_est = _usr_eso.vel_est  + v_est_dot*dt;
    _usr_eso.delta_est = _usr_eso.delta_est  + d_est_dot*dt;
}


void PosControl::getLocalPositionSetpoint(vehicle_local_position_setpoint_s &local_position_setpoint) const
{
        local_position_setpoint.x = _pos_sp(0);
        local_position_setpoint.y = _pos_sp(1);
        local_position_setpoint.z = _pos_sp(2);
        local_position_setpoint.yaw = _yaw_sp;
        local_position_setpoint.yawspeed = _yawspeed_sp;
        local_position_setpoint.vx = _vel_sp(0);
        local_position_setpoint.vy = _vel_sp(1);
        local_position_setpoint.vz = _vel_sp(2);
        _acc_sp.copyTo(local_position_setpoint.acceleration);
        _thr_sp.copyTo(local_position_setpoint.thrust);
}


void PosControl::getAttitudeSetpoint(vehicle_attitude_setpoint_s &attitude_setpoint) const
{
        ControlMath::thrustToAttitude(_thr_sp, _yaw_sp, attitude_setpoint);
        attitude_setpoint.yaw_sp_move_rate = _yawspeed_sp;
}

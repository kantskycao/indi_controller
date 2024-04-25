#include "Att_control.hpp"

using namespace matrix;

void addIfNotNan(float &setpoint, const float addition)
{
	if (PX4_ISFINITE(setpoint) && PX4_ISFINITE(addition)) {
		// No NAN, add to the setpoint
		setpoint += addition;

	} else if (!PX4_ISFINITE(setpoint)) {
		// Setpoint NAN, take addition
		setpoint = addition;
	}

	// Addition is NAN or both are NAN, nothing to do
}

void addIfNotNanVector3f(Vector3f &setpoint, const Vector3f &addition)
{
	for (int i = 0; i < 3; i++) {
		addIfNotNan(setpoint(i), addition(i));
	}
}

void setZeroIfNanVector3f(Vector3f &vector)
{
	// Adding zero vector overwrites elements that are NaN with zero
	addIfNotNanVector3f(vector, Vector3f());
}
void Att_Control::setMotorSpeed(const motor_speed_s &motor_speed)
{
    _motor_speed = motor_speed;
}

void Att_Control::setControllerGain(const matrix::Vector3f &KQ,const matrix::Vector3f &KW)
{
     _controller_param.K_xi(0,0)=KQ(0);
     _controller_param.K_xi(1,1)=KQ(1);
     _controller_param.K_xi(2,2)=KQ(2);

     _controller_param.K_omega(0,0)=KW(0);
     _controller_param.K_omega(1,1)=KW(1);
     _controller_param.K_omega(2,2)=KW(2);
}

void Att_Control::setpropParams(const float prop2force, const float prop2torque, const double arm_length){
    if(or_rotorcraft_max_rotors == 4){
        /*
        body frame
                X-axiz
                ^
        3       |      1
                |
        <-------0------>Y-axiz
                |
        2       |      4
        */
        _R_Torque(0,0) = -0.707 * arm_length * prop2force;
        _R_Torque(0,1) = 0.707 * arm_length * prop2force;
        _R_Torque(0,2) = 0.707 * arm_length * prop2force;
        _R_Torque(0,3) = -0.707 * arm_length * prop2force;
        _R_Torque(1,0) = 0.707 * arm_length * prop2force;
        _R_Torque(1,1) = -0.707 * arm_length * prop2force;
        _R_Torque(1,2) = 0.707 * arm_length * prop2force;
        _R_Torque(1,3) = -0.707 * arm_length * prop2force;
        _R_Torque(1,0) = 1 * prop2torque;
        _R_Torque(1,1) = 1 * prop2torque;
        _R_Torque(1,2) = -1 * prop2torque;
        _R_Torque(1,3) = -1 * prop2torque;
    }
    else if(or_rotorcraft_max_rotors == 6){
        /*
        body frame
            X-axiz
            ^
            3  |  5
            |
        2<-----0---->1  Y-axiz
            |
            6  |  4
        */
        _R_Torque(0,0) = -1 * arm_length * prop2force;
        _R_Torque(0,1) = 1 * arm_length * prop2force;
        _R_Torque(0,2) = 1 * 0.5 * arm_length * prop2force;
        _R_Torque(0,3) = -1 * 0.5 * arm_length * prop2force;
        _R_Torque(0,4) = -1 * 0.5 * arm_length * prop2force;
        _R_Torque(0,5) = 1 * 0.5 * arm_length * prop2force;
        _R_Torque(1,0) = 0;
        _R_Torque(1,1) = 0;
        _R_Torque(1,2) = 1 * 0.866 * arm_length * prop2force;
        _R_Torque(1,3) = -1 * 0.866 * arm_length * prop2force;
        _R_Torque(1,4) = 1 * 0.866 * arm_length * prop2force;
        _R_Torque(1,5) = -1 * 0.866 * arm_length * prop2force;
        _R_Torque(2,0) = -1 * prop2torque;
        _R_Torque(2,1) = 1 * prop2torque;
        _R_Torque(2,2) = -1 * prop2torque;
        _R_Torque(2,3) = 1 * prop2torque;
        _R_Torque(2,4) = 1 * prop2torque;
        _R_Torque(2,5) = -1 * prop2torque;
    }
    else
        _R_Torque.setZero();

}
void Att_Control:: update(const matrix::Quatf &q, const matrix::Vector3f &rate, const matrix::Vector3f &tau_static,
                          const float &dt, const bool &landed, Vector3f& torque,matrix::Vector3f &rates_sp)
{
    matrix::Vector3f  tau_static_cheack = tau_static;
    setZeroIfNanVector3f(tau_static_cheack);
    tau_static_cheack(0) = math::constrain(tau_static_cheack(0), -2.0f,2.0f);
    tau_static_cheack(1) = math::constrain(tau_static_cheack(1), -2.0f,2.0f);
    tau_static_cheack(2) = math::constrain(tau_static_cheack(2), -2.0f,2.0f);
    runAttitudeControl(q,rate,tau_static_cheack,dt,torque,rates_sp);
}

void  Att_Control::runAttitudeControl(const matrix::Quatf &q, const matrix::Vector3f &rate,
                                      const matrix::Vector3f &tau_static,const float &dt, Vector3f& torque,matrix::Vector3f &rates_sp)
{
    Eulerf euler(q);
    Eulerf euler_sp(_attitude_setpoint_q);

    double wprop2_arr[or_rotorcraft_max_rotors];
    for (int i = 0; i < or_rotorcraft_max_rotors; ++i) {
        wprop2_arr[i] = _motor_speed.motor_spd_rpm[i] * _motor_speed.motor_spd_rpm[i];
    }
    ema_filter_generic(wprop2_arr,or_rotorcraft_max_rotors, _ema_generic_filter.scale, or_rotorcraft_max_rotors,
                       or_rotorcraft_max_rotors, _ema_generic_filter.in, _ema_generic_filter.alpha, _ema_generic_filter.out);

    matrix::Matrix<float, or_rotorcraft_max_rotors, 1>  wprop2_vec;
    for (int i = 0; i < or_rotorcraft_max_rotors; ++i) {
        wprop2_vec.row(i) = _ema_generic_filter.out[i];
    }

    double wf_arr[3] = { rate(0), rate(1), rate(2) };   //这里使用的是已经滤波之后的角速度值，但也可以获取IMU的测量角速度，使用ema进行滤波。不知道使用哪一种更好
    ema_filter(wf_arr, _ema_filter.scale, _ema_filter.in, _ema_filter.alpha, _ema_filter.out);
    _wf = Vector3f(_ema_filter.out[0], _ema_filter.out[1], _ema_filter.out[2]);
    _dot_Omega_f = (_wf - _wf_prev) / dt; // 检查
    _wf_prev = _wf;

    //omega_ref 表示角速度参考值，dot_Omeag_ref 表示角加速度参考值，此处都设置为0
    Vector3f omega_ref = Vector3f(0.0f,0.0f,0.0f);
    Vector3f dot_Omeag_ref = Vector3f(0.0f,0.0f,0.0f);

    Vector3f dot_Omeag_c = _controller_param.K_xi * (euler_sp - euler) + _controller_param.K_omega * (omega_ref - rate) + dot_Omeag_ref;
    Vector3f mu_f = _R_Torque * wprop2_vec;
    torque = mu_f + _I_b * (dot_Omeag_c - _dot_Omega_f);// 和孙的代码不一样
    //##########
    _tau(0) = torque(0);
    _tau(1) = torque(1);
    _tau(2) = torque(2);
}

void Att_Control::UsrAttitudeESO(matrix::Vector3f bm_omega,matrix::Vector3f u,float dt)
{
    Matrix<float,3, 3> ESO_gain;
    ESO_gain.identity();
    ESO_gain(0,0)=_usr_eso.bm_gain(0);
    ESO_gain(1,1)=_usr_eso.bm_gain(1);
    ESO_gain(2,2)=_usr_eso.bm_gain(2);

    // check
    setZeroIfNanVector3f(_usr_eso.delta_esti);
    setZeroIfNanVector3f(_usr_eso.bm_rate_esti);
    setZeroIfNanVector3f(bm_omega);
    setZeroIfNanVector3f(u);
    dt = PX4_ISFINITE(dt)? dt:0.0f;

    _usr_eso.bm_rate_esti_dot=u + _usr_eso.delta_esti + 2.0f*ESO_gain*(bm_omega-_usr_eso.bm_rate_esti);
    _usr_eso.delta_esti_dot=ESO_gain*ESO_gain*(bm_omega-_usr_eso.bm_rate_esti);
    _usr_eso.bm_rate_esti=_usr_eso.bm_rate_esti+_usr_eso.bm_rate_esti_dot*dt;
    _usr_eso.delta_esti=_usr_eso.delta_esti+_usr_eso.delta_esti_dot*dt;
}

void Att_Control::resetESO()
{
    _usr_eso.bm_rate_esti={0.0f,0.0f,0.0f};
    _usr_eso.bm_rate_esti_dot={0.0f,0.0f,0.0f};
    _usr_eso.delta_esti={0.0f,0.0f,0.0f};
    _usr_eso.delta_esti_dot={0.0f,0.0f,0.0f};
    _tau={0.0f,0.0f,0.0f};
}

void Att_Control::getRateControlStatus(rate_ctrl_status_s &rate_ctrl_status)
{
    rate_ctrl_status.rollspeed_integ = _usr_eso.delta_esti(0);
    rate_ctrl_status.pitchspeed_integ = _usr_eso.delta_esti(1);
    rate_ctrl_status.yawspeed_integ = _usr_eso.delta_esti(2);
}

void Att_Control::setEMAFilterParams(){
    double scale_temp[9] = {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0};
    double in_temp[3] = { 0.0, 0.0, 0.0 };
    double alpha_temp[3] = { 0.001, 0.001, 0.001 };
    double out_temp[3] = { 0.0, 0.0, 0.0 };

    std::copy(scale_temp, scale_temp + 9, _ema_filter.scale);
    std::copy(in_temp, in_temp + 9, _ema_filter.in);
    std::copy(alpha_temp, alpha_temp + 9, _ema_filter.alpha);
    std::copy(out_temp, out_temp + 9, _ema_filter.out);
}

void Att_Control::setEMAFilterGenericParams(){
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
    std::copy(scale_temp, scale_temp + 9, _ema_generic_filter.scale);
    std::copy(in_temp, in_temp + 9, _ema_generic_filter.in);
    std::copy(alpha_temp, alpha_temp + 9, _ema_generic_filter.alpha);
    std::copy(out_temp, out_temp + 9, _ema_generic_filter.out);

}
/* --- Filter Definition ------------------------------------------------------- */
/* Exponential moving average filter(EMA) generic */
void Att_Control::ema_filter(const double raw[3], const double scale[3*3],
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

void Att_Control::ema_filter_generic(const double raw[], size_t raw_size,
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

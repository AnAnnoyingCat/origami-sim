#include "crease_constraints.h"
#include <iostream>

static Eigen::VectorXd edge_angle_prev;

void setup_prev_angle(int size){
	edge_angle_prev.resize(size);
	edge_angle_prev.setZero();
}

void F_crease(Eigen::Ref<Eigen::Matrix<double, 12, 1>> f, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> q4, double k_crease, double theta_target, int creaseID) {
	// Get all the relevant variables together
	double h1, h2;
	Eigen::Vector3d n1, n2;
	double alpha4_31, alpha3_14, alpha4_23, alpha3_42;
	
	geth(h1, q3, q4, q1);
	geth(h2, q3, q4, q2);
	getNormal(n1, q1, q4, q3);
	getNormal(n2, q2, q3, q4);
	getAngle(alpha4_31, q4, q3, q1);
	cot(alpha4_31);						// We only need the cot of the angles, so precompute it
	getAngle(alpha3_14, q3, q1, q4);
	cot(alpha3_14);
	getAngle(alpha4_23, q4, q2, q3);
	cot(alpha4_23);
	getAngle(alpha3_42, q3, q4, q2);
	cot(alpha3_42);

	// Calculate current fold angle theta 
    Eigen::Vector3d crease_dir = (q4 - q3).normalized();
	double current_theta = std::atan2((n1.cross(n2)).dot(crease_dir), n1.dot(n2));

	// Unwrap the angle
	double delta = current_theta - edge_angle_prev(creaseID);
	if (delta > M_PI){
		current_theta -= 2 * M_PI;
	} else if (delta < -M_PI){
		current_theta += 2 * M_PI;
	}
	edge_angle_prev(creaseID) = current_theta;

	// Precompute some values
	Eigen::Vector3d n1h1 = n1 / h1;
	Eigen::Vector3d n2h2 = n2 / h2;
	
	double kdtheta = -k_crease * (current_theta - theta_target);

	// dθ/dp1
	f.segment<3>(0) = kdtheta * n1 / h1;

	// dθ/dp2
	f.segment<3>(3) = kdtheta * n2 / h2;
	
	// dθ/dp3
	f.segment<3>(6) = kdtheta * ((-alpha4_31 / (alpha3_14 + alpha4_31)) * n1h1 + (-alpha4_23 / (alpha3_42 + alpha4_23)) * n2h2);

	// dθ/dp4
	f.segment<3>(9) = kdtheta * ((-alpha3_14 / (alpha3_14 + alpha4_31)) * n1h1 + (-alpha3_42 / (alpha3_42 + alpha4_23)) * n2h2);
}

void F_crease_maple(Eigen::Ref<Eigen::Matrix<double, 12, 1>> f, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> q4, double k_crease, double theta_target){
	f.setZero();

	// x[0] = q3; x[1] = q4, x[2] = q1, x[3] = q2

	// Handwritten
	double x[4][3];
	x[0][0] = q3(0);
	x[0][1] = q3(1);
	x[0][2] = q3(2);

	x[1][0] = q4(0);
	x[1][1] = q4(1);
	x[1][2] = q4(2);

	x[2][0] = q1(0);
	x[2][1] = q1(1);
	x[2][2] = q1(2);

	x[3][0] = q2(0);
	x[3][1] = q2(1);
	x[3][2] = q2(2);
	double fx[12];
	double MapleGenVar1, MapleGenVar2, MapleGenVar3, MapleGenVar4, MapleGenVar5;

	// Maple generated
	double t1 = -x[2][1]+x[1][1];
	double t2 = x[0][2]-x[2][2];
	double t4 = -x[2][2]+x[1][2];
	double t5 = x[0][1]-x[2][1];
	double t7 = t1*t2-t4*t5;
	double t8 = t7*t7;
	double t9 = x[1][0]-x[2][0];
	double t11 = x[0][0]-x[2][0];
	double t13 = t4*t11-t9*t2;
	double t14 = t13*t13;
	double t17 = -t1*t11+t9*t5;
	double t18 = t17*t17;
	double t19 = t8+t14+t18;
	double t20 = sqrt(t19);
	double t21 = 1/t20;
	double t22 = t21*t13;
	double t23 = x[0][1]-x[3][1];
	double t24 = x[1][2]-x[3][2];
	double t26 = x[0][2]-x[3][2];
	double t27 = x[1][1]-x[3][1];
	double t29 = t23*t24-t26*t27;
	double t30 = t29*t29;
	double t31 = x[0][0]-x[3][0];
	double t33 = x[1][0]-x[3][0];
	double t35 = -t31*t24+t26*t33;
	double t36 = t35*t35;
	double t39 = -t23*t33+t31*t27;
	double t40 = t39*t39;
	double t41 = t30+t36+t40;
	double t42 = sqrt(t41);
	double t43 = 1/t42;
	double t44 = t43*t39;
	double t47 = t21*t17;
	double t48 = t43*t35;
	double t51 = 0.1E1*t22*t44-0.1E1*t47*t48;
	double t52 = x[1][0]-x[0][0];
	double t53 = t52*t52;
	double t54 = x[1][1]-x[0][1];
	double t55 = t54*t54;
	double t56 = x[1][2]-x[0][2];
	double t57 = t56*t56;
	double t58 = t53+t55+t57;
	double t59 = sqrt(t58);
	double t60 = 1/t59;
	double t61 = t51*t60;
	double t64 = t21*t7;
	double t67 = t43*t29;
	double t70 = -0.1E1*t64*t44+0.1E1*t47*t67;
	double t71 = t70*t60;
	double t78 = 0.1E1*t64*t48-0.1E1*t22*t67;
	double t79 = t78*t60;
	double t82 = 0.1E1*t61*t52+0.1E1*t71*t54+0.1E1*t79*t56;
	double t89 = 0.1E1+0.1E1*t64*t67+0.1E1*t22*t48+0.1E1*t47*t44;
	double t90 = 1/t89;
	double t92 = atan(t82*t90);
	double t95 = k_crease*(2.0*t92-theta_target);
	double t97 = 1/t20/t19;
	double t98 = t97*t13;
	double t101 = -t17*t1+t13*t4;
	double t102 = 2.0*t44*t101;
	double t105 = t21*t4;
	double t109 = 1/t42/t41;
	double t110 = t109*t39;
	double t113 = -t35*t24+t39*t27;
	double t114 = 2.0*t110*t113;
	double t117 = t43*t27;
	double t120 = t97*t17;
	double t121 = 2.0*t48*t101;
	double t124 = -t21*t1;
	double t127 = t109*t35;
	double t128 = 2.0*t127*t113;
	double t131 = -t43*t24;
	double t139 = 1/t59/t58;
	double t140 = t51*t139;
	double t144 = 0.1E1*t61;
	double t145 = t97*t7;
	double t152 = 2.0*t67*t101;
	double t157 = t109*t29;
	double t158 = 2.0*t157*t113;
	double t165 = t70*t139;
	double t166 = -2.0*t54*t52;
	double t185 = t78*t139;
	double t186 = -2.0*t56*t52;
	double t191 = t89*t89;
	double t192 = 1/t191;
	double t193 = t82*t192;
	double t217 = t82*t82;
	double t220 = 1/(t217*t192+1.0);
	double t226 = t17*t9-t7*t4;
	double t227 = 2.0*t44*t226;
	double t232 = t29*t24-t39*t33;
	double t233 = 2.0*t110*t232;
	double t236 = -t43*t33;
	double t239 = 2.0*t48*t226;
	double t242 = t21*t9;
	double t245 = 2.0*t127*t232;
	double t256 = -t21*t4;
	double t263 = 2.0*t67*t226;
	double t268 = 2.0*t157*t232;
	double t271 = t43*t24;
	double t281 = 0.1E1*t71;
	double t298 = -2.0*t56*t54;
	double t331 = t7*t1-t13*t9;
	double t332 = 2.0*t44*t331;
	double t335 = -t21*t9;
	double t340 = -t29*t27+t35*t33;
	double t341 = 2.0*t110*t340;
	double t344 = 2.0*t48*t331;
	double t347 = 2.0*t127*t340;
	double t350 = t43*t33;
	double t361 = t21*t1;
	double t366 = 2.0*t67*t331;
	double t369 = 2.0*t157*t340;
	double t372 = -t43*t27;
	double t404 = 0.1E1*t79;
	double t435 = -t13*t2+t17*t5;
	double t436 = 2.0*t44*t435;
	double t439 = -t21*t2;
	double t444 = -t39*t23+t35*t26;
	double t445 = 2.0*t110*t444;
	double t448 = -t43*t23;
	double t451 = 2.0*t48*t435;
	double t454 = t21*t5;
	double t457 = 2.0*t127*t444;
	double t460 = t43*t26;
	double t476 = 2.0*t67*t435;
	double t481 = 2.0*t157*t444;
	double t488 = 2.0*t54*t52;
	double t507 = 2.0*t56*t52;
	double t540 = -t17*t11+t7*t2;
	double t541 = 2.0*t44*t540;
	double t546 = -t29*t26+t39*t31;
	double t547 = 2.0*t110*t546;
	double t550 = t43*t31;
	double t553 = 2.0*t48*t540;
	double t556 = -t21*t11;
	double t559 = 2.0*t127*t546;
	double t570 = t21*t2;
	double t577 = 2.0*t67*t540;
	double t582 = 2.0*t157*t546;
	double t585 = -t43*t26;
	double t611 = 2.0*t56*t54;
	double t644 = t13*t11-t7*t5;
	double t645 = 2.0*t44*t644;
	double t648 = t21*t11;
	double t653 = t29*t23-t35*t31;
	double t654 = 2.0*t110*t653;
	double t657 = 2.0*t48*t644;
	double t660 = 2.0*t127*t653;
	double t663 = -t43*t31;
	double t674 = -t21*t5;
	double t679 = 2.0*t67*t644;
	double t682 = 2.0*t157*t653;
	double t685 = t43*t23;
	double t747 = -t13*t56+t17*t54;
	double t748 = 2.0*t44*t747;
	double t751 = -t21*t56;
	double t754 = 2.0*t48*t747;
	double t757 = t21*t54;
	double t766 = 2.0*t67*t747;
	double t805 = -t17*t52+t7*t56;
	double t806 = 2.0*t44*t805;
	double t809 = 2.0*t48*t805;
	double t812 = -t21*t52;
	double t821 = t21*t56;
	double t824 = 2.0*t67*t805;
	double t863 = t13*t52-t7*t54;
	double t864 = 2.0*t44*t863;
	double t867 = t21*t52;
	double t870 = 2.0*t48*t863;
	double t879 = -t21*t54;
	double t882 = 2.0*t67*t863;
	double t921 = t35*t56-t39*t54;
	double t922 = 2.0*t110*t921;
	double t925 = -t43*t54;
	double t928 = 2.0*t127*t921;
	double t931 = t43*t56;
	double t942 = 2.0*t157*t921;
	double t979 = -t29*t56+t39*t52;
	double t980 = 2.0*t110*t979;
	double t983 = t43*t52;
	double t986 = 2.0*t127*t979;
	double t997 = 2.0*t157*t979;
	double t1000 = -t43*t56;
	double t1037 = t29*t54-t35*t52;
	double t1038 = 2.0*t110*t1037;
	double t1041 = 2.0*t127*t1037;
	double t1044 = -t43*t52;
	double t1053 = 2.0*t157*t1037;
	double t1056 = t43*t54;
	MapleGenVar1 = -0.2E1*t95;
	MapleGenVar4 = (0.1E1*(-0.5*t98*t102+0.1E1*t105*t44-0.5*t22*t114+0.1E1*
	t22*t117+0.5*t120*t121-0.1E1*t124*t48+0.5*t47*t128-0.1E1*t47*t131)*t60*t52+
	0.1E1*t140*t52*t52-t144+0.1E1*(0.5*t145*t102+0.5*t64*t114-0.1E1*t64*t117-0.5*
	t120*t152+0.1E1*t124*t67-0.5*t47*t158)*t60*t54-0.5*t165*t166+0.1E1*(-0.5*t145*
	t121-0.5*t64*t128+0.1E1*t64*t131+0.5*t98*t152-0.1E1*t105*t67+0.5*t22*t158)*t60*
	t56-0.5*t185*t186)*t90;
	MapleGenVar5 = -t193*(-0.5*t145*t152-0.5*t64*t158-0.5*t98*t121+0.1E1*t105
	*t48-0.5*t22*t128+0.1E1*t22*t131-0.5*t120*t102+0.1E1*t124*t44-0.5*t47*t114+
	0.1E1*t47*t117);
	MapleGenVar3 = MapleGenVar4+MapleGenVar5;
	MapleGenVar4 = t220;
	MapleGenVar2 = MapleGenVar3*MapleGenVar4;
	fx[0] = MapleGenVar1*MapleGenVar2;
	MapleGenVar1 = -0.2E1*t95;
	MapleGenVar4 = (0.1E1*(-0.5*t98*t227-0.5*t22*t233+0.1E1*t22*t236+0.5*t120
	*t239-0.1E1*t242*t48+0.5*t47*t245)*t60*t52-0.5*t140*t166+0.1E1*(0.5*t145*t227
	-0.1E1*t256*t44+0.5*t64*t233-0.1E1*t64*t236-0.5*t120*t263+0.1E1*t242*t67-0.5*
	t47*t268+0.1E1*t47*t271)*t60*t54+0.1E1*t165*t54*t54-t281+0.1E1*(-0.5*t145*t239+
	0.1E1*t256*t48-0.5*t64*t245+0.5*t98*t263+0.5*t22*t268-0.1E1*t22*t271)*t60*t56
	-0.5*t185*t298)*t90;
	MapleGenVar5 = -t193*(-0.5*t145*t263+0.1E1*t256*t67-0.5*t64*t268+0.1E1*
	t64*t271-0.5*t98*t239-0.5*t22*t245-0.5*t120*t227+0.1E1*t242*t44-0.5*t47*t233+
	0.1E1*t47*t236);
	MapleGenVar3 = MapleGenVar4+MapleGenVar5;
	MapleGenVar4 = t220;
	MapleGenVar2 = MapleGenVar3*MapleGenVar4;
	fx[1] = MapleGenVar1*MapleGenVar2;
	MapleGenVar1 = -0.2E1*t95;
	MapleGenVar4 = (0.1E1*(-0.5*t98*t332+0.1E1*t335*t44-0.5*t22*t341+0.5*t120
	*t344+0.5*t47*t347-0.1E1*t47*t350)*t60*t52-0.5*t140*t186+0.1E1*(0.5*t145*t332
	-0.1E1*t361*t44+0.5*t64*t341-0.5*t120*t366-0.5*t47*t369+0.1E1*t47*t372)*t60*t54
	-0.5*t165*t298+0.1E1*(-0.5*t145*t344+0.1E1*t361*t48-0.5*t64*t347+0.1E1*t64*t350
	+0.5*t98*t366-0.1E1*t335*t67+0.5*t22*t369-0.1E1*t22*t372)*t60*t56+0.1E1*t185*
	t56*t56-t404)*t90;
	MapleGenVar5 = -t193*(-0.5*t145*t366+0.1E1*t361*t67-0.5*t64*t369+0.1E1*
	t64*t372-0.5*t98*t344+0.1E1*t335*t48-0.5*t22*t347+0.1E1*t22*t350-0.5*t120*t332
	-0.5*t47*t341);
	MapleGenVar3 = MapleGenVar4+MapleGenVar5;
	MapleGenVar4 = t220;
	MapleGenVar2 = MapleGenVar3*MapleGenVar4;
	fx[2] = MapleGenVar1*MapleGenVar2;
	MapleGenVar1 = -0.2E1*t95;
	MapleGenVar4 = (0.1E1*(-0.5*t98*t436+0.1E1*t439*t44-0.5*t22*t445+0.1E1*
	t22*t448+0.5*t120*t451-0.1E1*t454*t48+0.5*t47*t457-0.1E1*t47*t460)*t60*t52
	-0.1E1*t140*t52*t52+t144+0.1E1*(0.5*t145*t436+0.5*t64*t445-0.1E1*t64*t448-0.5*
	t120*t476+0.1E1*t454*t67-0.5*t47*t481)*t60*t54-0.5*t165*t488+0.1E1*(-0.5*t145*
	t451-0.5*t64*t457+0.1E1*t64*t460+0.5*t98*t476-0.1E1*t439*t67+0.5*t22*t481)*t60*
	t56-0.5*t185*t507)*t90;
	MapleGenVar5 = -t193*(-0.5*t145*t476-0.5*t64*t481-0.5*t98*t451+0.1E1*t439
	*t48-0.5*t22*t457+0.1E1*t22*t460-0.5*t120*t436+0.1E1*t454*t44-0.5*t47*t445+
	0.1E1*t47*t448);
	MapleGenVar3 = MapleGenVar4+MapleGenVar5;
	MapleGenVar4 = t220;
	MapleGenVar2 = MapleGenVar3*MapleGenVar4;
	fx[3] = MapleGenVar1*MapleGenVar2;
	MapleGenVar1 = -0.2E1*t95;
	MapleGenVar4 = (0.1E1*(-0.5*t98*t541-0.5*t22*t547+0.1E1*t22*t550+0.5*t120
	*t553-0.1E1*t556*t48+0.5*t47*t559)*t60*t52-0.5*t140*t488+0.1E1*(0.5*t145*t541
	-0.1E1*t570*t44+0.5*t64*t547-0.1E1*t64*t550-0.5*t120*t577+0.1E1*t556*t67-0.5*
	t47*t582+0.1E1*t47*t585)*t60*t54-0.1E1*t165*t54*t54+t281+0.1E1*(-0.5*t145*t553+
	0.1E1*t570*t48-0.5*t64*t559+0.5*t98*t577+0.5*t22*t582-0.1E1*t22*t585)*t60*t56
	-0.5*t185*t611)*t90;
	MapleGenVar5 = -t193*(-0.5*t145*t577+0.1E1*t570*t67-0.5*t64*t582+0.1E1*
	t64*t585-0.5*t98*t553-0.5*t22*t559-0.5*t120*t541+0.1E1*t556*t44-0.5*t47*t547+
	0.1E1*t47*t550);
	MapleGenVar3 = MapleGenVar4+MapleGenVar5;
	MapleGenVar4 = t220;
	MapleGenVar2 = MapleGenVar3*MapleGenVar4;
	fx[4] = MapleGenVar1*MapleGenVar2;
	MapleGenVar1 = -0.2E1*t95;
	MapleGenVar4 = (0.1E1*(-0.5*t98*t645+0.1E1*t648*t44-0.5*t22*t654+0.5*t120
	*t657+0.5*t47*t660-0.1E1*t47*t663)*t60*t52-0.5*t140*t507+0.1E1*(0.5*t145*t645
	-0.1E1*t674*t44+0.5*t64*t654-0.5*t120*t679-0.5*t47*t682+0.1E1*t47*t685)*t60*t54
	-0.5*t165*t611+0.1E1*(-0.5*t145*t657+0.1E1*t674*t48-0.5*t64*t660+0.1E1*t64*t663
	+0.5*t98*t679-0.1E1*t648*t67+0.5*t22*t682-0.1E1*t22*t685)*t60*t56-0.1E1*t185*
	t56*t56+t404)*t90;
	MapleGenVar5 = -t193*(-0.5*t145*t679+0.1E1*t674*t67-0.5*t64*t682+0.1E1*
	t64*t685-0.5*t98*t657+0.1E1*t648*t48-0.5*t22*t660+0.1E1*t22*t663-0.5*t120*t645
	-0.5*t47*t654);
	MapleGenVar3 = MapleGenVar4+MapleGenVar5;
	MapleGenVar4 = t220;
	MapleGenVar2 = MapleGenVar3*MapleGenVar4;
	fx[5] = MapleGenVar1*MapleGenVar2;
	fx[6] = -0.2E1*t95*((0.1E1*(-0.5*t98*t748+0.1E1*t751*t44+0.5*t120*t754
	-0.1E1*t757*t48)*t60*t52+0.1E1*(0.5*t145*t748-0.5*t120*t766+0.1E1*t757*t67)*t60
	*t54+0.1E1*(-0.5*t145*t754+0.5*t98*t766-0.1E1*t751*t67)*t60*t56)*t90-t193*(-0.5
	*t145*t766-0.5*t98*t754+0.1E1*t751*t48-0.5*t120*t748+0.1E1*t757*t44))*t220;
		fx[7] = -0.2E1*t95*((0.1E1*(-0.5*t98*t806+0.5*t120*t809-0.1E1*t812*t48)*
	t60*t52+0.1E1*(0.5*t145*t806-0.1E1*t821*t44-0.5*t120*t824+0.1E1*t812*t67)*t60*
	t54+0.1E1*(-0.5*t145*t809+0.1E1*t821*t48+0.5*t98*t824)*t60*t56)*t90-t193*(-0.5*
	t145*t824+0.1E1*t821*t67-0.5*t98*t809-0.5*t120*t806+0.1E1*t812*t44))*t220;
		fx[8] = -0.2E1*t95*((0.1E1*(-0.5*t98*t864+0.1E1*t867*t44+0.5*t120*t870)*
	t60*t52+0.1E1*(0.5*t145*t864-0.1E1*t879*t44-0.5*t120*t882)*t60*t54+0.1E1*(-0.5*
	t145*t870+0.1E1*t879*t48+0.5*t98*t882-0.1E1*t867*t67)*t60*t56)*t90-t193*(-0.5*
	t145*t882+0.1E1*t879*t67-0.5*t98*t870+0.1E1*t867*t48-0.5*t120*t864))*t220;
		fx[9] = -0.2E1*t95*((0.1E1*(-0.5*t22*t922+0.1E1*t22*t925+0.5*t47*t928
	-0.1E1*t47*t931)*t60*t52+0.1E1*(0.5*t64*t922-0.1E1*t64*t925-0.5*t47*t942)*t60*
	t54+0.1E1*(-0.5*t64*t928+0.1E1*t64*t931+0.5*t22*t942)*t60*t56)*t90-t193*(-0.5*
	t64*t942-0.5*t22*t928+0.1E1*t22*t931-0.5*t47*t922+0.1E1*t47*t925))*t220;
		fx[10] = -0.2E1*t95*((0.1E1*(-0.5*t22*t980+0.1E1*t22*t983+0.5*t47*t986)*
	t60*t52+0.1E1*(0.5*t64*t980-0.1E1*t64*t983-0.5*t47*t997+0.1E1*t47*t1000)*t60*
	t54+0.1E1*(-0.5*t64*t986+0.5*t22*t997-0.1E1*t22*t1000)*t60*t56)*t90-t193*(-0.5*
	t64*t997+0.1E1*t64*t1000-0.5*t22*t986-0.5*t47*t980+0.1E1*t47*t983))*t220;
		fx[11] = -0.2E1*t95*((0.1E1*(-0.5*t22*t1038+0.5*t47*t1041-0.1E1*t47*t1044
	)*t60*t52+0.1E1*(0.5*t64*t1038-0.5*t47*t1053+0.1E1*t47*t1056)*t60*t54+0.1E1*(
	-0.5*t64*t1041+0.1E1*t64*t1044+0.5*t22*t1053-0.1E1*t22*t1056)*t60*t56)*t90-t193
	*(-0.5*t64*t1053+0.1E1*t64*t1056-0.5*t22*t1041+0.1E1*t22*t1044-0.5*t47*t1038))*
	t220;

	// x[0] = q3; x[1] = q4, x[2] = q1, x[3] = q2
	f(0) = fx[6];
	f(1) = fx[7];
	f(2) = fx[8];
	f(3) = fx[9];
	f(4) = fx[10];
	f(5) = fx[11];
	f(6) = fx[0];
	f(7) = fx[1];
	f(8) = fx[2];
	f(9) = fx[3];
	f(10) = fx[4];
	f(11) = fx[5];

	for (int i = 0; i < 12; i++){
		f(i) = fx[i];
	}

	f *= -1;
}
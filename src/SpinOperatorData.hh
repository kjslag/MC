// developer: Kevin Slagle (kslagle@physics.ucsb.edu)
// date: 2013

SpinOperator_<1> ising_flipper({ [](Spin_<1> &s){s = -s;} });

const Float sqrt1_2 = sqrt(.5 );
const Float sqrt3_4 = sqrt(.75);
#define OP(...) [](Spin_<COUNT_ARGS(__VA_ARGS__)> &s) { s = {{__VA_ARGS__}}; }
#define r sqrt1_2
#define q sqrt3_4

SpinOperator_<2> signed_permutation_flipper_2({ OP(-s[0],s[1]), OP(-s[1],-s[0]), OP(s[1],s[0]), OP(s[0],-s[1]) });

SpinOperator_<3> signed_permutation_flipper_3({
OP(-s[0], s[1], s[2]), OP(-s[1],-s[0], s[2]), OP(-s[2], s[1],-s[0]),
OP( s[2], s[1], s[0]), OP( s[1], s[0], s[2]), OP( s[0],-s[1], s[2]),
OP( s[0],-s[2],-s[1]), OP( s[0], s[2], s[1]), OP( s[0], s[1],-s[2]) });

SpinOperator_<4> signed_permutation_flipper_4({
OP(-s[0], s[1], s[2], s[3]), OP(-s[1],-s[0], s[2], s[3]),
OP(-s[2], s[1],-s[0], s[3]), OP(-s[3], s[1], s[2],-s[0]),
OP( s[3], s[1], s[2], s[0]), OP( s[2], s[1], s[0], s[3]),
OP( s[1], s[0], s[2], s[3]), OP( s[0],-s[1], s[2], s[3]),
OP( s[0],-s[2],-s[1], s[3]), OP( s[0],-s[3], s[2],-s[1]),
OP( s[0], s[3], s[2], s[1]), OP( s[0], s[2], s[1], s[3]),
OP( s[0], s[1],-s[2], s[3]), OP( s[0], s[1],-s[3],-s[2]),
OP( s[0], s[1], s[3], s[2]), OP( s[0], s[1], s[2],-s[3]) });

SpinOperator_<2> cos2_flipper({ OP(-s[0],s[1]), OP( s[0],-s[1]) });
SpinOperator_<2> cos3_flipper({ OP(s[0],-s[1]), OP(-s[0]/2-q*s[1],-q*s[0]+s[1]/2), OP(-s[0]/2+q*s[1],q*s[0]+s[1]/2) }, false);
SpinOperator_<2> cos4_flipper({ OP(-s[0],s[1]), OP(-s[1],-s[0]), OP(s[1],s[0]), OP(s[0],-s[1]) });

SpinOperator_<2> cos6_flipper({
OP(-s[0], s[1]),  OP( s[0],-s[1]),
OP(-s[0]/2-q*s[1],-q*s[0]+s[1]/2), OP(-s[0]/2+q*s[1], q*s[0]+s[1]/2),
OP( s[0]/2-q*s[1],-q*s[0]-s[1]/2), OP( s[0]/2+q*s[1], q*s[0]-s[1]/2) });

SpinOperator_<2> cos8_flipper({
OP(-s[0], s[1]),  OP(-s[1],-s[0]), OP( s[1], s[0]),  OP( s[0],-s[1]),
OP(-r*s[0]-r*s[1],-r*s[0]+r*s[1]), OP(-r*s[0]+r*s[1], r*s[0]+r*s[1]),
OP( r*s[0]-r*s[1],-r*s[0]-r*s[1]), OP( r*s[0]+r*s[1], r*s[0]-r*s[1]) });

SpinOperator_<6> vison_triangle_flipper({OP(-s[0],-s[1], s[2], s[3], s[4], s[5]),
OP(-s[0], s[1], s[2], s[3],-s[5],-s[4]),
OP(-s[0], s[1], s[2], s[3], s[5], s[4]),
OP(-s[1],-s[0],-s[2], s[3], s[4], s[5]),
OP(-s[1],-s[0], s[2],-s[3], s[4], s[5]),
OP( s[1], s[0],-s[2], s[3], s[4], s[5]),
OP( s[1], s[0], s[2],-s[3], s[4], s[5]),
OP( s[0],-s[1], s[2], s[3],-s[5],-s[4]),
OP( s[0],-s[1], s[2], s[3], s[5], s[4]),
OP( s[0], s[1],-s[2],-s[3], s[4], s[5]),
OP( s[0], s[1],-s[3],-s[2],-s[4], s[5]),
OP( s[0], s[1],-s[3],-s[2], s[4],-s[5]),
OP( s[0], s[1], s[3], s[2],-s[4], s[5]),
OP( s[0], s[1], s[3], s[2], s[4],-s[5]),
OP( s[0], s[1], s[2], s[3],-s[4],-s[5]),
OP(-r*s[4]-r*s[5],-r*s[4]+r*s[5],-r*s[2]-r*s[3],-r*s[2]+r*s[3],-r*s[0]-r*s[1],-r*s[0]+r*s[1]),
OP(-r*s[4]-r*s[5],-r*s[4]+r*s[5], r*s[2]+r*s[3], r*s[2]-r*s[3],-r*s[0]-r*s[1],-r*s[0]+r*s[1]),
OP(-r*s[4]-r*s[5], r*s[4]-r*s[5],-r*s[2]-r*s[3],-r*s[2]+r*s[3],-r*s[0]+r*s[1],-r*s[0]-r*s[1]),
OP(-r*s[4]-r*s[5], r*s[4]-r*s[5], r*s[2]+r*s[3], r*s[2]-r*s[3],-r*s[0]+r*s[1],-r*s[0]-r*s[1]),
OP(-r*s[4]+r*s[5],-r*s[4]-r*s[5],-r*s[2]+r*s[3], r*s[2]+r*s[3],-r*s[0]-r*s[1], r*s[0]-r*s[1]),
OP(-r*s[4]+r*s[5],-r*s[4]-r*s[5], r*s[2]-r*s[3],-r*s[2]-r*s[3],-r*s[0]-r*s[1], r*s[0]-r*s[1]),
OP(-r*s[4]+r*s[5], r*s[4]+r*s[5],-r*s[2]+r*s[3], r*s[2]+r*s[3],-r*s[0]+r*s[1], r*s[0]+r*s[1]),
OP(-r*s[4]+r*s[5], r*s[4]+r*s[5], r*s[2]-r*s[3],-r*s[2]-r*s[3],-r*s[0]+r*s[1], r*s[0]+r*s[1]),
OP( r*s[4]-r*s[5],-r*s[4]-r*s[5],-r*s[2]+r*s[3], r*s[2]+r*s[3], r*s[0]-r*s[1],-r*s[0]-r*s[1]),
OP( r*s[4]-r*s[5],-r*s[4]-r*s[5], r*s[2]-r*s[3],-r*s[2]-r*s[3], r*s[0]-r*s[1],-r*s[0]-r*s[1]),
OP( r*s[4]-r*s[5], r*s[4]+r*s[5],-r*s[2]+r*s[3], r*s[2]+r*s[3], r*s[0]+r*s[1],-r*s[0]+r*s[1]),
OP( r*s[4]-r*s[5], r*s[4]+r*s[5], r*s[2]-r*s[3],-r*s[2]-r*s[3], r*s[0]+r*s[1],-r*s[0]+r*s[1]),
OP( r*s[4]+r*s[5],-r*s[4]+r*s[5],-r*s[2]-r*s[3],-r*s[2]+r*s[3], r*s[0]-r*s[1], r*s[0]+r*s[1]),
OP( r*s[4]+r*s[5],-r*s[4]+r*s[5], r*s[2]+r*s[3], r*s[2]-r*s[3], r*s[0]-r*s[1], r*s[0]+r*s[1]),
OP( r*s[4]+r*s[5], r*s[4]-r*s[5],-r*s[2]-r*s[3],-r*s[2]+r*s[3], r*s[0]+r*s[1], r*s[0]-r*s[1]),
OP( r*s[4]+r*s[5], r*s[4]-r*s[5], r*s[2]+r*s[3], r*s[2]-r*s[3], r*s[0]+r*s[1], r*s[0]-r*s[1]),
OP(-r*s[2]-r*s[3],-r*s[2]+r*s[3],-r*s[0]-r*s[1],-r*s[0]+r*s[1],-r*s[4]-r*s[5],-r*s[4]+r*s[5]),
OP(-r*s[2]-r*s[3],-r*s[2]+r*s[3],-r*s[0]-r*s[1],-r*s[0]+r*s[1], r*s[4]+r*s[5], r*s[4]-r*s[5]),
OP(-r*s[2]-r*s[3], r*s[2]-r*s[3],-r*s[0]+r*s[1],-r*s[0]-r*s[1],-r*s[4]+r*s[5], r*s[4]+r*s[5]),
OP(-r*s[2]-r*s[3], r*s[2]-r*s[3],-r*s[0]+r*s[1],-r*s[0]-r*s[1], r*s[4]-r*s[5],-r*s[4]-r*s[5]),
OP(-r*s[2]+r*s[3],-r*s[2]-r*s[3],-r*s[0]-r*s[1], r*s[0]-r*s[1],-r*s[4]-r*s[5],-r*s[4]+r*s[5]),
OP(-r*s[2]+r*s[3],-r*s[2]-r*s[3],-r*s[0]-r*s[1], r*s[0]-r*s[1], r*s[4]+r*s[5], r*s[4]-r*s[5]),
OP(-r*s[2]+r*s[3], r*s[2]+r*s[3],-r*s[0]+r*s[1], r*s[0]+r*s[1],-r*s[4]+r*s[5], r*s[4]+r*s[5]),
OP(-r*s[2]+r*s[3], r*s[2]+r*s[3],-r*s[0]+r*s[1], r*s[0]+r*s[1], r*s[4]-r*s[5],-r*s[4]-r*s[5]),
OP( r*s[2]-r*s[3],-r*s[2]-r*s[3], r*s[0]-r*s[1],-r*s[0]-r*s[1],-r*s[4]+r*s[5], r*s[4]+r*s[5]),
OP( r*s[2]-r*s[3],-r*s[2]-r*s[3], r*s[0]-r*s[1],-r*s[0]-r*s[1], r*s[4]-r*s[5],-r*s[4]-r*s[5]),
OP( r*s[2]-r*s[3], r*s[2]+r*s[3], r*s[0]+r*s[1],-r*s[0]+r*s[1],-r*s[4]-r*s[5],-r*s[4]+r*s[5]),
OP( r*s[2]-r*s[3], r*s[2]+r*s[3], r*s[0]+r*s[1],-r*s[0]+r*s[1], r*s[4]+r*s[5], r*s[4]-r*s[5]),
OP( r*s[2]+r*s[3],-r*s[2]+r*s[3], r*s[0]-r*s[1], r*s[0]+r*s[1],-r*s[4]+r*s[5], r*s[4]+r*s[5]),
OP( r*s[2]+r*s[3],-r*s[2]+r*s[3], r*s[0]-r*s[1], r*s[0]+r*s[1], r*s[4]-r*s[5],-r*s[4]-r*s[5]),
OP( r*s[2]+r*s[3], r*s[2]-r*s[3], r*s[0]+r*s[1], r*s[0]-r*s[1],-r*s[4]-r*s[5],-r*s[4]+r*s[5]),
OP( r*s[2]+r*s[3], r*s[2]-r*s[3], r*s[0]+r*s[1], r*s[0]-r*s[1], r*s[4]+r*s[5], r*s[4]-r*s[5]),
OP(-r*s[0]-r*s[1],-r*s[0]+r*s[1],-r*s[4]-r*s[5],-r*s[4]+r*s[5],-r*s[2]-r*s[3],-r*s[2]+r*s[3]),
OP(-r*s[0]-r*s[1],-r*s[0]+r*s[1],-r*s[4]+r*s[5],-r*s[4]-r*s[5],-r*s[2]-r*s[3], r*s[2]-r*s[3]),
OP(-r*s[0]-r*s[1],-r*s[0]+r*s[1], r*s[4]-r*s[5], r*s[4]+r*s[5], r*s[2]+r*s[3],-r*s[2]+r*s[3]),
OP(-r*s[0]-r*s[1],-r*s[0]+r*s[1], r*s[4]+r*s[5], r*s[4]-r*s[5], r*s[2]+r*s[3], r*s[2]-r*s[3]),
OP(-r*s[0]+r*s[1], r*s[0]+r*s[1],-r*s[4]-r*s[5], r*s[4]-r*s[5],-r*s[2]+r*s[3],-r*s[2]-r*s[3]),
OP(-r*s[0]+r*s[1], r*s[0]+r*s[1],-r*s[4]+r*s[5], r*s[4]+r*s[5],-r*s[2]+r*s[3], r*s[2]+r*s[3]),
OP(-r*s[0]+r*s[1], r*s[0]+r*s[1], r*s[4]-r*s[5],-r*s[4]-r*s[5], r*s[2]-r*s[3],-r*s[2]-r*s[3]),
OP(-r*s[0]+r*s[1], r*s[0]+r*s[1], r*s[4]+r*s[5],-r*s[4]+r*s[5], r*s[2]-r*s[3], r*s[2]+r*s[3]),
OP( r*s[0]-r*s[1],-r*s[0]-r*s[1],-r*s[4]-r*s[5], r*s[4]-r*s[5],-r*s[2]+r*s[3],-r*s[2]-r*s[3]),
OP( r*s[0]-r*s[1],-r*s[0]-r*s[1],-r*s[4]+r*s[5], r*s[4]+r*s[5],-r*s[2]+r*s[3], r*s[2]+r*s[3]),
OP( r*s[0]-r*s[1],-r*s[0]-r*s[1], r*s[4]-r*s[5],-r*s[4]-r*s[5], r*s[2]-r*s[3],-r*s[2]-r*s[3]),
OP( r*s[0]-r*s[1],-r*s[0]-r*s[1], r*s[4]+r*s[5],-r*s[4]+r*s[5], r*s[2]-r*s[3], r*s[2]+r*s[3]),
OP( r*s[0]+r*s[1], r*s[0]-r*s[1],-r*s[4]-r*s[5],-r*s[4]+r*s[5],-r*s[2]-r*s[3],-r*s[2]+r*s[3]),
OP( r*s[0]+r*s[1], r*s[0]-r*s[1],-r*s[4]+r*s[5],-r*s[4]-r*s[5],-r*s[2]-r*s[3], r*s[2]-r*s[3]),
OP( r*s[0]+r*s[1], r*s[0]-r*s[1], r*s[4]-r*s[5], r*s[4]+r*s[5], r*s[2]+r*s[3],-r*s[2]+r*s[3]),
OP( r*s[0]+r*s[1], r*s[0]-r*s[1], r*s[4]+r*s[5], r*s[4]-r*s[5], r*s[2]+r*s[3], r*s[2]-r*s[3]) });

#undef OP
#undef r
#undef q

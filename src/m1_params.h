void get_m1params_v3(float cushion, float tmax, float *moments, float *m1_params, double gmax, bool debug)
{
    int i_cushion = (int) (cushion*100 + .01);
    int i_tmax = (int) (tmax*100 + .01);
    int i_gmax = (int) (gmax + .01);

    std::cout << i_cushion << " -- " << i_tmax << " -- " << i_gmax << std::endl;

    if ((i_cushion == 75) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -0.768176;
        m1_params[1] = 3.487883;
        m1_params[2] = -0.453589;
        m1_params[3] = 2.401463;
        m1_params[4] = 1.364426;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 80) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -0.790123;
        m1_params[1] = 3.604938;
        m1_params[2] = -0.395062;
        m1_params[3] = 2.255144;
        m1_params[4] = 1.514403;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 85) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -0.592593;
        m1_params[1] = 3.407407;
        m1_params[2] = -0.444444;
        m1_params[3] = 2.222222;
        m1_params[4] = 1.777778;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 90) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -0.888889;
        m1_params[1] = 4.000000;
        m1_params[2] = 0.000000;
        m1_params[3] = 2.222222;
        m1_params[4] = 1.777778;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 75) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -0.592593;
        m1_params[1] = 3.473251;
        m1_params[2] = -0.746228;
        m1_params[3] = 2.335620;
        m1_params[4] = 1.759488;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 80) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = 0.153635;
        m1_params[1] = 3.429355;
        m1_params[2] = -1.163237;
        m1_params[3] = 2.189300;
        m1_params[4] = 2.139918;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 85) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -0.592593;
        m1_params[1] = 3.407407;
        m1_params[2] = -0.444444;
        m1_params[3] = 2.222222;
        m1_params[4] = 1.777778;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 90) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -0.888889;
        m1_params[1] = 4.000000;
        m1_params[2] = 0.000000;
        m1_params[3] = 2.222222;
        m1_params[4] = 1.777778;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 75) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = -0.310115;
        m1_params[1] = 3.687040;
        m1_params[2] = -1.028299;
        m1_params[3] = 2.274247;
        m1_params[4] = 1.970025;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 80) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = -0.708022;
        m1_params[1] = 3.587055;
        m1_params[2] = -0.474724;
        m1_params[3] = 2.336432;
        m1_params[4] = 1.677996;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 85) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = -0.438957;
        m1_params[1] = 3.517147;
        m1_params[2] = -0.724280;
        m1_params[3] = 2.233196;
        m1_params[4] = 1.986283;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 90) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = 0.000000;
        m1_params[1] = 3.604938;
        m1_params[2] = -0.790123;
        m1_params[3] = 2.255144;
        m1_params[4] = 2.008230;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 75) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -0.518214;
        m1_params[1] = 3.576487;
        m1_params[2] = -1.054311;
        m1_params[3] = 2.482345;
        m1_params[4] = 2.113093;
        m1_params[5] = 32.000000;
    } else if ((i_cushion == 80) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -0.708022;
        m1_params[1] = 3.587055;
        m1_params[2] = -0.474724;
        m1_params[3] = 2.336432;
        m1_params[4] = 1.677996;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 85) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = 0.131687;
        m1_params[1] = 3.604938;
        m1_params[2] = -1.207133;
        m1_params[3] = 2.145405;
        m1_params[4] = 2.161866;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 90) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -0.592593;
        m1_params[1] = 3.407407;
        m1_params[2] = -0.592593;
        m1_params[3] = 2.123457;
        m1_params[4] = 1.679012;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 75) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -0.016664;
        m1_params[1] = 3.879693;
        m1_params[2] = -1.361988;
        m1_params[3] = 2.364883;
        m1_params[4] = 2.571356;
        m1_params[5] = 33.000000;
    } else if ((i_cushion == 80) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -0.438957;
        m1_params[1] = 3.648834;
        m1_params[2] = -0.537723;
        m1_params[3] = 2.540466;
        m1_params[4] = 2.013717;
        m1_params[5] = 33.000000;
    } else if ((i_cushion == 85) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -0.234111;
        m1_params[1] = 3.502515;
        m1_params[2] = -1.123000;
        m1_params[3] = 2.496571;
        m1_params[4] = 2.138089;
        m1_params[5] = 32.000000;
    } else if ((i_cushion == 90) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -0.221917;
        m1_params[1] = 3.500889;
        m1_params[2] = -0.855154;
        m1_params[3] = 2.244577;
        m1_params[4] = 2.019408;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 75) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -0.185744;
        m1_params[1] = 4.113804;
        m1_params[2] = -1.082355;
        m1_params[3] = 2.312859;
        m1_params[4] = 2.649393;
        m1_params[5] = 35.000000;
    } else if ((i_cushion == 80) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -0.438957;
        m1_params[1] = 3.648834;
        m1_params[2] = -0.537723;
        m1_params[3] = 2.540466;
        m1_params[4] = 2.013717;
        m1_params[5] = 33.000000;
    } else if ((i_cushion == 85) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -0.234111;
        m1_params[1] = 3.502515;
        m1_params[2] = -0.966926;
        m1_params[3] = 2.158411;
        m1_params[4] = 1.982015;
        m1_params[5] = 32.000000;
    } else if ((i_cushion == 90) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -0.221917;
        m1_params[1] = 3.448865;
        m1_params[2] = -0.803130;
        m1_params[3] = 2.244577;
        m1_params[4] = 2.019408;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 75) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = 0.065031;
        m1_params[1] = 3.388711;
        m1_params[2] = -0.039018;
        m1_params[3] = 2.840014;
        m1_params[4] = 2.077122;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 80) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = 0.206066;
        m1_params[1] = 3.947975;
        m1_params[2] = -0.964894;
        m1_params[3] = 2.976172;
        m1_params[4] = 2.850988;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 85) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = 0.514556;
        m1_params[1] = 3.778895;
        m1_params[2] = -1.162018;
        m1_params[3] = 2.522583;
        m1_params[4] = 2.655083;
        m1_params[5] = 35.000000;
    } else if ((i_cushion == 90) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = 0.318651;
        m1_params[1] = 3.749632;
        m1_params[2] = -1.337195;
        m1_params[3] = 1.981202;
        m1_params[4] = 2.440482;
        m1_params[5] = 35.000000;
    } else if ((i_cushion == 75) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = 0.065031;
        m1_params[1] = 3.388711;
        m1_params[2] = 0.013006;
        m1_params[3] = 2.840014;
        m1_params[4] = 2.077122;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 80) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = 0.102017;
        m1_params[1] = 3.843926;
        m1_params[2] = -1.120967;
        m1_params[3] = 3.015191;
        m1_params[4] = 2.850988;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 85) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = 0.159325;
        m1_params[1] = 3.840675;
        m1_params[2] = -1.048214;
        m1_params[3] = 2.468933;
        m1_params[4] = 2.534776;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 90) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = 0.409694;
        m1_params[1] = 3.619570;
        m1_params[2] = -0.868973;
        m1_params[3] = 2.228319;
        m1_params[4] = 2.534776;
        m1_params[5] = 34.000000;
    }

}















void get_m1params_v3_cardiac22(float cushion, float tmax, float *moments, float *m1_params, double gmax, bool debug)
{
    int i_cushion = (int) (cushion*100 + .01);
    int i_tmax = (int) (tmax*100 + .01);
    int i_gmax = (int) (gmax + .01);

    std::cout << i_cushion << " -- " << i_tmax << " -- " << i_gmax << std::endl;

    if ((i_cushion == 75) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -0.490169;
        m1_params[1] = 3.261088;
        m1_params[2] = -0.709648;
        m1_params[3] = 2.339278;
        m1_params[4] = 1.752172;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 80) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = 0.000000;
        m1_params[1] = 2.814815;
        m1_params[2] = -0.790123;
        m1_params[3] = 2.518519;
        m1_params[4] = 1.876543;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 85) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -0.296296;
        m1_params[1] = 3.012346;
        m1_params[2] = -0.395062;
        m1_params[3] = 2.222222;
        m1_params[4] = 1.185185;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 90) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -0.888889;
        m1_params[1] = 3.111111;
        m1_params[2] = 0.000000;
        m1_params[3] = 2.222222;
        m1_params[4] = 0.888889;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 75) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -0.526749;
        m1_params[1] = 3.209877;
        m1_params[2] = -0.592593;
        m1_params[3] = 2.178326;
        m1_params[4] = 1.569273;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 80) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -0.263374;
        m1_params[1] = 3.209877;
        m1_params[2] = -0.790123;
        m1_params[3] = 2.123457;
        m1_params[4] = 1.481481;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 85) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -1.283951;
        m1_params[1] = 3.506173;
        m1_params[2] = -0.395062;
        m1_params[3] = 2.024691;
        m1_params[4] = 1.580247;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 90) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -0.888889;
        m1_params[1] = 3.407407;
        m1_params[2] = -0.296296;
        m1_params[3] = 2.222222;
        m1_params[4] = 1.777778;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 75) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = 0.078037;
        m1_params[1] = 3.316364;
        m1_params[2] = -1.068130;
        m1_params[3] = 2.232383;
        m1_params[4] = 2.398212;
        m1_params[5] = 32.000000;
    } else if ((i_cushion == 80) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = 0.060966;
        m1_params[1] = 3.250521;
        m1_params[2] = -0.975461;
        m1_params[3] = 2.400650;
        m1_params[4] = 2.247015;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 85) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = 0.000000;
        m1_params[1] = 3.209877;
        m1_params[2] = -0.921811;
        m1_params[3] = 2.320988;
        m1_params[4] = 1.876543;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 90) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = 0.592593;
        m1_params[1] = 3.078189;
        m1_params[2] = -1.119342;
        m1_params[3] = 2.320988;
        m1_params[4] = 2.271605;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 75) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = 0.042270;
        m1_params[1] = 3.211502;
        m1_params[2] = -1.202662;
        m1_params[3] = 2.277092;
        m1_params[4] = 2.303104;
        m1_params[5] = 32.000000;
    } else if ((i_cushion == 80) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = 0.008942;
        m1_params[1] = 3.276533;
        m1_params[2] = -1.105523;
        m1_params[3] = 2.348626;
        m1_params[4] = 2.194991;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 85) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = 0.219479;
        m1_params[1] = 3.209877;
        m1_params[2] = -1.207133;
        m1_params[3] = 2.101509;
        m1_params[4] = 2.227709;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 90) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = 0.592593;
        m1_params[1] = 3.078189;
        m1_params[2] = -1.119342;
        m1_params[3] = 2.320988;
        m1_params[4] = 2.271605;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 75) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -0.171112;
        m1_params[1] = 3.639079;
        m1_params[2] = -1.208352;
        m1_params[3] = 2.267337;
        m1_params[4] = 2.194178;
        m1_params[5] = 35.000000;
    } else if ((i_cushion == 80) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = 0.774679;
        m1_params[1] = 3.393588;
        m1_params[2] = -1.680638;
        m1_params[3] = 2.150282;
        m1_params[4] = 2.945283;
        m1_params[5] = 33.000000;
    } else if ((i_cushion == 85) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = 0.125997;
        m1_params[1] = 3.404969;
        m1_params[2] = -0.742570;
        m1_params[3] = 2.646141;
        m1_params[4] = 2.101509;
        m1_params[5] = 32.000000;
    } else if ((i_cushion == 90) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -0.175583;
        m1_params[1] = 3.341564;
        m1_params[2] = -0.888889;
        m1_params[3] = 2.160037;
        m1_params[4] = 2.081390;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 75) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = 0.466189;
        m1_params[1] = 3.463496;
        m1_params[2] = -1.559518;
        m1_params[3] = 2.049484;
        m1_params[4] = 2.558350;
        m1_params[5] = 35.000000;
    } else if ((i_cushion == 80) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -0.301580;
        m1_params[1] = 3.365950;
        m1_params[2] = -0.913275;
        m1_params[3] = 2.293350;
        m1_params[4] = 2.164914;
        m1_params[5] = 33.000000;
    } else if ((i_cushion == 85) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = 0.125997;
        m1_params[1] = 3.404969;
        m1_params[2] = -0.742570;
        m1_params[3] = 2.646141;
        m1_params[4] = 2.101509;
        m1_params[5] = 32.000000;
    } else if ((i_cushion == 90) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -0.091043;
        m1_params[1] = 3.170858;
        m1_params[2] = -0.847025;
        m1_params[3] = 2.290098;
        m1_params[4] = 2.035869;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 75) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = 0.254433;
        m1_params[1] = 3.511457;
        m1_params[2] = -0.145506;
        m1_params[3] = 3.004217;
        m1_params[4] = 2.178123;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 80) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = -0.162577;
        m1_params[1] = 3.078189;
        m1_params[2] = -0.315399;
        m1_params[3] = 3.242798;
        m1_params[4] = 2.085251;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 85) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = -0.344663;
        m1_params[1] = 2.468526;
        m1_params[2] = -0.266626;
        m1_params[3] = 2.274247;
        m1_params[4] = 2.375451;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 90) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = -0.367424;
        m1_params[1] = 2.883097;
        m1_params[2] = -0.263374;
        m1_params[3] = 2.716049;
        m1_params[4] = 1.991770;
        m1_params[5] = 35.000000;
    } else if ((i_cushion == 75) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = 0.254433;
        m1_params[1] = 3.511457;
        m1_params[2] = -0.145506;
        m1_params[3] = 3.004217;
        m1_params[4] = 2.178123;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 80) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = -0.162577;
        m1_params[1] = 3.078189;
        m1_params[2] = -0.315399;
        m1_params[3] = 3.242798;
        m1_params[4] = 2.085251;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 85) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = -0.344663;
        m1_params[1] = 2.780674;
        m1_params[2] = -0.266626;
        m1_params[3] = 2.248235;
        m1_params[4] = 2.375451;
        m1_params[5] = 35.000000;
    } else if ((i_cushion == 90) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = -0.367424;
        m1_params[1] = 2.883097;
        m1_params[2] = -0.263374;
        m1_params[3] = 2.716049;
        m1_params[4] = 1.991770;
        m1_params[5] = 35.000000;
    }

}







void get_m1params_v3_neuro(float cushion, float tmax, float *moments, float *m1_params, double gmax, bool debug)
{
    int i_cushion = (int) (cushion*100 + .01);
    int i_tmax = (int) (tmax*100 + .01);
    int i_gmax = (int) (gmax + .01);

    std::cout << i_cushion << " -- " << i_tmax << " -- " << i_gmax << std::endl;

    if ((i_cushion == 75) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -3.661840;
        m1_params[1] = 8.308896;
        m1_params[2] = 0.026012;
        m1_params[3] = 4.390184;
        m1_params[4] = 1.793019;
        m1_params[5] = 40.000000;
    } else if ((i_cushion == 80) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -3.947975;
        m1_params[1] = 8.000000;
        m1_params[2] = 0.104049;
        m1_params[3] = 4.026012;
        m1_params[4] = 2.258396;
        m1_params[5] = 39.000000;
    } else if ((i_cushion == 85) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -3.954478;
        m1_params[1] = 8.130061;
        m1_params[2] = -0.149571;
        m1_params[3] = 4.250368;
        m1_params[4] = 1.872479;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 90) && (i_tmax == 80) && (i_gmax == 45)) {
        m1_params[0] = -3.863435;
        m1_params[1] = 7.895951;
        m1_params[2] = 0.227608;
        m1_params[3] = 4.042270;
        m1_params[4] = 1.401006;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 75) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -4.130061;
        m1_params[1] = 8.556013;
        m1_params[2] = -0.260123;
        m1_params[3] = 4.260123;
        m1_params[4] = 2.092161;
        m1_params[5] = 39.000000;
    } else if ((i_cushion == 80) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -3.863435;
        m1_params[1] = 8.000000;
        m1_params[2] = -0.104049;
        m1_params[3] = 3.895951;
        m1_params[4] = 1.964131;
        m1_params[5] = 39.000000;
    } else if ((i_cushion == 85) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -3.512269;
        m1_params[1] = 7.986994;
        m1_params[2] = -0.409694;
        m1_params[3] = 4.045522;
        m1_params[4] = 1.820454;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 90) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -4.006503;
        m1_params[1] = 7.739877;
        m1_params[2] = 0.227608;
        m1_params[3] = 4.676320;
        m1_params[4] = 1.453031;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 75) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = -3.313926;
        m1_params[1] = 8.939694;
        m1_params[2] = 0.052025;
        m1_params[3] = 5.090484;
        m1_params[4] = 2.449830;
        m1_params[5] = 42.000000;
    } else if ((i_cushion == 80) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = -3.904080;
        m1_params[1] = 8.544632;
        m1_params[2] = 0.299141;
        m1_params[3] = 4.752731;
        m1_params[4] = 2.252299;
        m1_params[5] = 41.000000;
    } else if ((i_cushion == 85) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = -3.492760;
        m1_params[1] = 8.149571;
        m1_params[2] = -0.123558;
        m1_params[3] = 4.253620;
        m1_params[4] = 2.468120;
        m1_params[5] = 39.000000;
    } else if ((i_cushion == 90) && (i_tmax == 70) && (i_gmax == 45)) {
        m1_params[0] = -4.000000;
        m1_params[1] = 8.000000;
        m1_params[2] = 0.156074;
        m1_params[3] = 4.052025;
        m1_params[4] = 1.333333;
        m1_params[5] = 39.000000;
    } else if ((i_cushion == 75) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -4.318651;
        m1_params[1] = 8.809226;
        m1_params[2] = 0.803942;
        m1_params[3] = 5.153889;
        m1_params[4] = 1.780420;
        m1_params[5] = 44.000000;
    } else if ((i_cushion == 80) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -3.895951;
        m1_params[1] = 8.000000;
        m1_params[2] = 0.000000;
        m1_params[3] = 4.000000;
        m1_params[4] = 3.407407;
        m1_params[5] = 40.000000;
    } else if ((i_cushion == 85) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -3.869939;
        m1_params[1] = 7.954478;
        m1_params[2] = 0.058528;
        m1_params[3] = 3.850429;
        m1_params[4] = 1.804806;
        m1_params[5] = 40.000000;
    } else if ((i_cushion == 90) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -4.000000;
        m1_params[1] = 8.052025;
        m1_params[2] = 0.156074;
        m1_params[3] = 4.052025;
        m1_params[4] = 1.333333;
        m1_params[5] = 39.000000;
    } else if ((i_cushion == 75) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -1.622720;
        m1_params[1] = 8.361733;
        m1_params[2] = 0.199766;
        m1_params[3] = 4.467408;
        m1_params[4] = 3.827669;
        m1_params[5] = 47.000000;
    } else if ((i_cushion == 80) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -2.302698;
        m1_params[1] = 8.579383;
        m1_params[2] = -1.185185;
        m1_params[3] = 4.221105;
        m1_params[4] = 3.837423;
        m1_params[5] = 45.000000;
    } else if ((i_cushion == 85) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -3.349489;
        m1_params[1] = 8.938272;
        m1_params[2] = 0.338160;
        m1_params[3] = 5.281309;
        m1_params[4] = 2.562617;
        m1_params[5] = 44.000000;
    } else if ((i_cushion == 90) && (i_tmax == 60) && (i_gmax == 45)) {
        m1_params[0] = -4.299141;
        m1_params[1] = 8.755982;
        m1_params[2] = 0.000000;
        m1_params[3] = 4.455215;
        m1_params[4] = 1.333333;
        m1_params[5] = 44.000000;
    } else if ((i_cushion == 75) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -0.256059;
        m1_params[1] = 9.296957;
        m1_params[2] = -1.839760;
        m1_params[3] = 4.130061;
        m1_params[4] = 4.555810;
        m1_params[5] = 44.000000;
    } else if ((i_cushion == 80) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -2.405934;
        m1_params[1] = 8.447086;
        m1_params[2] = -0.041457;
        m1_params[3] = 4.521872;
        m1_params[4] = 2.958086;
        m1_params[5] = 44.000000;
    } else if ((i_cushion == 85) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -3.056851;
        m1_params[1] = 9.196769;
        m1_params[2] = 0.351166;
        m1_params[3] = 5.099223;
        m1_params[4] = 3.251943;
        m1_params[5] = 44.000000;
    } else if ((i_cushion == 90) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -3.074938;
        m1_params[1] = 8.771427;
        m1_params[2] = -0.334908;
        m1_params[3] = 4.792969;
        m1_params[4] = 2.846924;
        m1_params[5] = 41.000000;
    } else if ((i_cushion == 75) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = -0.273129;
        m1_params[1] = 7.713865;
        m1_params[2] = -0.117055;
        m1_params[3] = 6.549611;
        m1_params[4] = 7.707362;
        m1_params[5] = 50.000000;
    } else if ((i_cushion == 80) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = -0.062592;
        m1_params[1] = 10.313062;
        m1_params[2] = -3.020475;
        m1_params[3] = 6.532541;
        m1_params[4] = 6.374841;
        m1_params[5] = 52.000000;
    } else if ((i_cushion == 85) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = -0.175583;
        m1_params[1] = 9.715592;
        m1_params[2] = -1.618656;
        m1_params[3] = 4.175583;
        m1_params[4] = 5.470508;
        m1_params[5] = 45.000000;
    } else if ((i_cushion == 90) && (i_tmax == 50) && (i_gmax == 45)) {
        m1_params[0] = -1.489407;
        m1_params[1] = 8.728547;
        m1_params[2] = -1.403445;
        m1_params[3] = 4.172535;
        m1_params[4] = 5.907230;
        m1_params[5] = 48.000000;
    } else if ((i_cushion == 75) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = -0.390184;
        m1_params[1] = 8.284103;
        m1_params[2] = -0.174770;
        m1_params[3] = 5.321343;
        m1_params[4] = 6.396992;
        m1_params[5] = 50.000000;
    } else if ((i_cushion == 80) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = 0.108926;
        m1_params[1] = 9.668445;
        m1_params[2] = -1.886704;
        m1_params[3] = 5.291267;
        m1_params[4] = 4.983996;
        m1_params[5] = 51.000000;
    } else if ((i_cushion == 85) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = -0.175583;
        m1_params[1] = 9.715592;
        m1_params[2] = -1.618656;
        m1_params[3] = 4.175583;
        m1_params[4] = 5.470508;
        m1_params[5] = 45.000000;
    } else if ((i_cushion == 90) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = -1.190469;
        m1_params[1] = 8.878321;
        m1_params[2] = 0.543210;
        m1_params[3] = 6.209623;
        m1_params[4] = 3.422649;
        m1_params[5] = 47.000000;
    }

}










void get_m1params_v3_ph160(float cushion, float tmax, float *moments, float *m1_params, double gmax, bool debug)
{
    int i_cushion = (int) (cushion*100 + .01);
    int i_tmax = (int) (tmax*100 + .01);
    int i_gmax = (int) (gmax + .01);

    std::cout << i_cushion << " -- " << i_tmax << " -- " << i_gmax << std::endl;

    if ((i_cushion == 90) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = 0.000000;
        m1_params[1] = 2.518519;
        m1_params[2] = -0.296296;
        m1_params[3] = 2.222222;
        m1_params[4] = 1.777778;
        m1_params[5] = 29.000000;
    } else if ((i_cushion == 90) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -0.592593;
        m1_params[1] = 1.925926;
        m1_params[2] = -0.395062;
        m1_params[3] = 2.716049;
        m1_params[4] = 1.481481;
        m1_params[5] = 30.000000;
    } else if ((i_cushion == 90) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = 0.065031;
        m1_params[1] = 2.397399;
        m1_params[2] = -0.764111;
        m1_params[3] = 2.385206;
        m1_params[4] = 1.940761;
        m1_params[5] = 31.000000;
    } else if ((i_cushion == 90) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = 0.618605;
        m1_params[1] = 2.024285;
        m1_params[2] = -0.836864;
        m1_params[3] = 2.232383;
        m1_params[4] = 2.115328;
        m1_params[5] = 34.000000;
    }

}







void get_m1params_v3_ph80(float cushion, float tmax, float *moments, float *m1_params, double gmax, bool debug)
{
    int i_cushion = (int) (cushion*100 + .01);
    int i_tmax = (int) (tmax*100 + .01);
    int i_gmax = (int) (gmax + .01);

    std::cout << i_cushion << " -- " << i_tmax << " -- " << i_gmax << std::endl;

    if ((i_cushion == 90) && (i_tmax == 80) && (i_gmax == 80)) {
        m1_params[0] = -0.556013;
        m1_params[1] = 2.662399;
        m1_params[2] = -0.156074;
        m1_params[3] = 2.510593;
        m1_params[4] = 1.852157;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 90) && (i_tmax == 70) && (i_gmax == 80)) {
        m1_params[0] = -0.951075;
        m1_params[1] = 2.602652;
        m1_params[2] = -0.058528;
        m1_params[3] = 2.613626;
        m1_params[4] = 1.774120;
        m1_params[5] = 38.000000;
    } else if ((i_cushion == 90) && (i_tmax == 60) && (i_gmax == 80)) {
        m1_params[0] = -0.416197;
        m1_params[1] = 2.783925;
        m1_params[2] = -0.111365;
        m1_params[3] = 2.764416;
        m1_params[4] = 2.419347;
        m1_params[5] = 39.000000;
    } else if ((i_cushion == 90) && (i_tmax == 50) && (i_gmax == 80)) {
        m1_params[0] = 1.691002;
        m1_params[1] = 2.537621;
        m1_params[2] = -0.573287;
        m1_params[3] = 4.274755;
        m1_params[4] = 3.748616;
        m1_params[5] = 44.000000;
    }

}
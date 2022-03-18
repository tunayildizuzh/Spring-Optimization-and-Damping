#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> 
using namespace std;

const double pi = 3.14159;
const double g = 9.8;

const double mgw = 0; // Weight of the wheel.
const double mgs = 0; // Weight of the spring.



// Position function of the spring for a lightly damped system.
double position_function(double A, double b, double m, double k, double t) { // Lightly damped.
    /* 
        F = m * a
        m*a + k*x + b * dx/dt = 0
        m * d2x/dt^2 + k*x + b*dx.dt = 0
        mx'' + bx' + kx = 0
        w = sqrt(k/m)
    */
   double a = k/m;
   double c = pow(b,2) / 4 * pow(m,2);
   double w_prime = sqrt(a - c);
   return A * exp((-b*t)/2*m) * cos(w_prime);
}

// Path in terms of a function. 
double path(double x) {
    //return sin(pi*x) * cos(x);
    return sin(0.4*x) * 5*cos(x);
}

// Retrieves y values of the path function. 
vector<double> y_val(double step_max, double increment) {
    vector<double> y_val;
    double h = 0.0;
    while(h < step_max) {
        y_val.push_back(path(h));
        h +=increment;
    }
    return y_val;
}

// Takes two point of y value and divides it with x step_size so we can get the angle at the end.
vector<double> angle(vector<double> y_val, double step_size) {
    double diff;
    vector<double> angle;
    for(int i = 0; i < y_val.size(); ++i) {
        diff = (y_val[i+1] - y_val[i]) / step_size;
        angle.push_back(atan(diff) * 180/pi);
    }
    return angle; 
}

// This checks when the angle returns from + to - or - to + which gives the exact vertex. 
bool sign_check(double num1, double num2) {
    bool status = false;

    if(num1 > 0) {
        if(num2 < 0) {
            status = true; 
        }
    }

    if(num1 < 0) {
        if(num2 > 0) {
            status = true;
        }
    }
    return status;
}


// Gets the vertex.
vector<double> vertex(vector<double> y_val) {
    vector<double> angles = angle(y_val, 0.001);
    vector<double> vertex;
    for(int i = 0; i < angles.size(); i++) {
        if(abs(angles[i]) < 1.0 ) {
            if(sign_check(angles[i],angles[i+1])) {
            vertex.push_back(y_val[i]);
            }
        }
    }
    return vertex;
}

double vertex_max(vector<double> vertex) { 
    vector<double> vertexes = vertex;
    double vertex_max;
    vertexes[1] = vertex_max;
    for(int i = 0; i < vertexes.size(); i++) {
        if(vertexes[i] > vertex_max) {
            vertex_max = vertexes[i];
        }
    }
    return vertex_max;
}

double vertex_min(vector<double> vertex) { 
    vector<double> vertexes = vertex;
    double vertex_min;
    vertexes[1] = vertex_min;
    for(int i = 0; i < vertexes.size(); i++) {
        if(vertexes[i] < vertex_min) {
            vertex_min = vertexes[i];
        }
    }
    return vertex_min;
}

vector<double> disp(vector<double> y_vals,double vertex_max) {
    double disp_max = vertex_max + 0.5; 
    double disp_neutral = (disp_max + vertex_min(vertex(y_val(20,0.001)))) / 2;
    vector<double> disp; 
    for(int i =0; i < y_vals.size(); ++i) {
        if(y_vals[i] < disp_neutral) { 
            disp.push_back(y_vals[i] - disp_neutral);
        } else if(y_vals[i] > disp_neutral) {
            disp.push_back(y_vals[i] - disp_neutral);
        }
    }
    return disp;
}

vector<double> spring_force(double kspring) {
    vector<double> sforce;
    vector<double> disps = disp(y_val(20,0.001), vertex_max(vertex(y_val(20,0.001))));
    for(int i = 0; i < disps.size(); ++i) {
        sforce.push_back(disps[i] * kspring);
    }
    return sforce;
}

vector<double> first_derivative(double ts_max, double h) {
    vector<double> f_der;
    for(int i =0; i*h < ts_max; ++i) {
        f_der.push_back((path(i*h + h) - path(i*h))/h);
    }

    return f_der;
}

vector<double> second_derivative(double ts_max, double h) {
    vector<double> s_der;
    // double path(double x) {
    // return sin(pi*x) * cos(x);
    // }
    for(int i =0; i*h < ts_max; ++i) {
       s_der.push_back((path(i*h + h) + path(i*h - h) - 2*path(i*h))/ pow(h,2));
    } 
    return s_der;
}

vector<double> radius(vector<double> f_der, vector<double> s_der) {
    vector<double> radius;
    double R;
    for(int i = 0; i < f_der.size(); ++i) {
        R = (pow((1+pow(f_der[i],2)),1.5))/(abs(s_der[i]));
        radius.push_back(R);
    }

    return radius;
}

vector<double> normal_force(vector<double> y_val, double m_spring, double m_wheel, vector<double> f_spring, double velocity, vector<double> radius) {
    vector<double> normal_force;
    vector<double> angles = angle(y_val, 0.001);
    vector<double> s_ders = second_derivative(20,0.001);
    vector<double> f_friction;
    f_friction.push_back(0.6 * 54.0837);
    //double f_friction [sizeof(s_ders)/sizeof(s_ders[0]) +1] = {0.6*54.0837};
    //f_friction[0] = 0.6 * 54.0837; // First normal force value without friction force.
    
    for(int i =0; i < s_ders.size(); ++i) {

        if(s_ders[i] > 0) {
            if(angles[i] > 0) {
                //normal_force.push_back((m_spring + m_wheel + f_spring[i] + ((m_spring+m_wheel) * pow(velocity,2)/radius[i]))/ cos(abs(angles[i]) * pi/180));
                normal_force.push_back( (((m_spring * g) + (m_wheel * g) + f_spring[i]) * cos(angles[i] * (pi/180)) + (f_friction[i] * sin(angles[i] * (pi/180)))) + ( ((m_spring+m_wheel) * pow(velocity,2)) / radius[i]) );
            }
            else {
                normal_force.push_back( (((m_spring * g) + (m_wheel * g) + f_spring[i]) * cos(angles[i] * (pi/180)) - (f_friction[i] * sin(angles[i] * (pi/180)))) + ( ((m_spring+m_wheel) * pow(velocity,2)) / radius[i]) );
            }
        }

        else {
            if(angles[i] > 0) {
                //normal_force.push_back((m_spring + m_wheel + f_spring[i] - ((m_spring+m_wheel) * pow(velocity,2)/radius[i]))/ cos(abs(angles[i]) * pi/180));
                normal_force.push_back( (((m_spring * g) + (m_wheel * g) + f_spring[i]) * cos(angles[i] * (pi/180)) + (f_friction[i] * sin(angles[i] * (pi/180)))) - ( ((m_spring+m_wheel) * pow(velocity,2)) / radius[i]) );
            }
            else {
                normal_force.push_back( (((m_spring * g) + (m_wheel * g) + f_spring[i]) * cos(angles[i] * (pi/180)) - (f_friction[i] * sin(angles[i] * (pi/180)))) - ( ((m_spring+m_wheel) * pow(velocity,2)) / radius[i]) );

            }
        }
        f_friction.push_back(0.6 * normal_force[i]);
    }
    
    return normal_force;
}

vector<double> c_occurs;
vector<double> c_loses;
bool contact_points(double num1, double num2, double num3) {
    bool status = false;
    if( num1 > 0) {
        if(num2 < 0) {
            c_loses.push_back(num3);
            //cout << "At x: " << num3 << " Contact loses. " << endl;
            //cout << "Contact loses at these points: " << endl;
            //cout << num3 << endl;
           
            status = true;
        }
    }

    if(num1 < 0) {
        if( num2 > 0) {
            c_occurs.push_back(num3);
            //cout << "At x: " << num3 << " Contact occurs. " << endl;
            cout << num3 << endl;
            status = true;
        }
    } 

    return status;
}

void display(int num) {
    vector<double> rads = radius(first_derivative(20,0.001),second_derivative(20,0.001));
    vector<double> n_forces = normal_force(y_val(20,0.001), 5, 15, spring_force(50), 2.7, rads);
    vector<double> s_derv = second_derivative(20,0.001);
    vector<double> f_spring = spring_force(40);
    vector<double> angles = angle(y_val(20,0.001), 0.001);
    vector<double> disps = disp(y_val(20,0.001), vertex_max(vertex(y_val(20,0.001))));
    switch(num) {
        case 1:   // NORMAL FORCE. 
            for(int i = 0; i < n_forces.size(); ++i) {
                cout << "At x: " << i * 0.001 << " Normal Force is: " << n_forces[i] << endl; 
            }
            break;
    
        case 2:  // CONTACT POINTS.
            for(int i=0; i < n_forces.size(); ++i) {
                if(contact_points(n_forces[i],n_forces[i+1], i*0.001)) {}
            }
            break;

        case 3:  // SPRING FORCE.
            for(int i = 0; i < f_spring.size(); ++i) {
                cout << "At x: " << i * 0.001 << " Spring Force is: " << f_spring[i] << endl;
            }
            break;

        case 4:  // SPRING DISPLACEMENT.
            for(int i = 0; i < disps.size(); ++i) {
                cout << "At x: " << i * 0.001 << " Displacement is: " << disps[i] << endl; 
            }
            break;
           
        case 5:  // ANGLES.
            for(int i = 0; i < angles.size(); ++i) {
                cout << "At x: " << i * 0.001 << "Angle is: " << angles[i] << endl;
            }
            break;
        case 6: 
            for(int i =0; i < n_forces.size(); ++i) {
                if(contact_points(n_forces[i], n_forces[i+1], i*0.001)) {
                } 
            }
            // for(double val: c_occurs) {
            //     cout << val << endl;
            // }
            for(double val: c_loses) {
                cout << val << endl;
            }
        case 7: // SPRING DISP RIGHT BEFORE TYRE LOSES CONTACT.
            for(int i =0; i < c_loses.size(); ++i) {
                cout<< "At x = " << c_loses[i] << "Spring displacement is: " << disps[c_loses[i]] << endl;

            }

    }


}


// (dx/dt)'' + 2 * C * w0 * (dx/dt) + w0 * x = 0
// C = c / 2* sqrt(m*k) 
// w0 = sqrt(k/m)

/* To use Runge Kutta, the ODE above can be written in terms of: 
    (dx/dt) = z
    (dx/dt)'' = dz/dt

    so the final form simplifies to 

    (dz/dt) + 2*C*w0*z + w0*x = 0;

    float c, m, k, C, w0, z, x;

    w0 = sqrt(k/m);
    C = c/(2*(sqrt(m*k)));
    
    return - 2 * C * w0 * z + w0 * x; 
*/


float f_func(float z) {
    return z;

};

float g_func(float x, float z) {
    float c, m, k, C, w0;
    // float c;
    // float m;
    // float k;
    // float C;
    // float w0;
    // dy/dt = (t,y)
    // 
    k = 50;
    m = 20;
    c = 200;

    w0 = sqrt(k/m);
    C = c / (2*(sqrt(m*k))); //Damping ratio.
    float result_new = -((k/m)*x) - ((c/m) * z);
    return result_new; 

}

float runga_kutta(float x0, float z0, const float h, float iter_max) {
    vector<float> rung;
    float x, z, k1, k2, k3, k4, l1,l2,l3,l4, k_total, l_total;
    float increment = 0;
    //float iter_max = 20;

    while(increment < iter_max) {
        k1 = h * f_func(z0);
        l1 = h * g_func(x0, z0);
        //cout << k1 << " k1 and l1 " << l1 << endl;

        k2 = h * f_func(z0 + (0.5*l1));
        l2 = h * g_func(x0 + (0.5*k1), (z0+(0.5*l1)));
        //cout << k2 << " k2 and l2 " << l2 << endl;

        k3 = h * f_func(z0 + (0.5*l2));
        l3 = h * g_func(x0 + (0.5*k2), (z0+(0.5*l2)));
        //cout << k3 << " k3 and l3 " << l3 << endl;

        k4 = h * f_func(z0 + l3);
        l4 = h * g_func((x0 + k3), (z0 + l3));
        //cout << k4 << " k4 and l4 " << l4 << endl;

        k_total = k1 + (2*k2) + (2*k3) + k4;
        l_total = l1 + (2*l2) + (2*k3) + l4;
        //cout << "k_tot " << k_total << endl;
        //cout << "l_tot " << l_total << endl;

        x = (x0 + ((0.16666)*(k_total)));
        //cout << x << endl;
        z = (z0 + ((0.16666)*(l_total)));
        //cout << z << endl;
        rung.push_back(x);
        x0 = x;
        z0 = z;
        increment += h;
        
    }

    for(float val: rung) {
        cout << val << endl;
    }
    return x;
}




int main() {
    vector<double> rads = radius(first_derivative(20,0.001),second_derivative(20,0.001));
    vector<double> n_forces = normal_force(y_val(20,0.001), 5, 15, spring_force(50), 2.7, rads);
    vector<double> s_derv = second_derivative(20,0.001);
    vector<double> f_spring = spring_force(40);
    vector<double> angles = angle(y_val(20,0.001), 0.001);
    vector<double> disps = disp(y_val(20,0.001), vertex_max(vertex(y_val(20,0.001))));
    vector<double> y_vals = y_val(20,0.001);
    vector<double> disps_contact;
    vector<double> flight_time;

    for(int i =0; i < n_forces.size(); ++i) {
            if(contact_points(n_forces[i], n_forces[i+1], i*0.001)) {
            } 
    }
    
    for(int i =0; i < c_loses.size(); ++i) {
        //cout<< "At x = " << c_loses[i] << "Spring displacement is: " << disps[c_loses[i]] << endl;
        disps_contact.push_back(disps[c_loses[i]]);
    }
    for(int i =0; i<c_loses.size(); ++i) {
        flight_time.push_back(c_occurs[i] - c_loses[i]);

    }
    flight_time.pop_back();
    for(float val:flight_time) {
        //cout << val << endl;
    }

    for(int i=0; i<flight_time.size(); ++i) {
        cout << i << " th flight," << "Damping datas as shown: " << endl;
        runga_kutta(disps_contact[i],2,0.01,12);
        cout << "END OF " << i << " th flight." << endl;
    }

    //display(2);
    //runga_kutta(1,0,0.01,12);
    // for(double val: y_vals) {
    //     cout << val << endl;
    // }

    // fstream pracdata2;
    // pracdata2.open("pracdata2.csv", ios::out);
    // if(pracdata2.is_open()) {
    //     pracdata2 << "Force values" << endl;
    //     for(double val: y_vals) {
    //     pracdata2 << val << endl;
    //     }
    //     pracdata2.close();
    // }
    
}

// energy. Visuals. damping. if i wanted to better what i would do.
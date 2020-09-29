#ifndef KMUCS_GRAPHICS_TRANSFORM_HPP
#define KMUCS_GRAPHICS_TRANSFORM_HPP

#include <cmath>
#include "vec.hpp"
#include "mat.hpp"
#include "operator.hpp"

namespace kmuvcl
{
    namespace math
    {
#ifndef M_PI
        const float M_PI = 3.14159265358979323846f;
#endif

        template <typename T>
        mat<4, 4, T> translate(T dx, T dy, T dz)
        {
            mat<4, 4, T> translateMat;

            // TODO: Fill up this function properly 
            for(int i = 0; i < 4; i++){
                translateMat[5*i] = 1;
            }
            translateMat[12] = dx;
            translateMat[13] = dy;
            translateMat[14] = dz;

            return translateMat;
        }

        template <typename T>
        mat<4, 4, T> rotate(T angle, T x, T y, T z)
        {
            mat<4, 4, T> rotateMat;

            // TODO: Fill up this function properly 
            double radian = (M_PI/180)*angle;
            double len, ux, uy, uz, cos_calc;
            len = sqrt(x*x+y*y+z*z);
            ux = x/len;
            uy = y/len;
            uz = z/len;
            cos_calc = 1-cos(radian);

            rotateMat[0] = cos(radian) + ux*ux*cos_calc;
            rotateMat[1] = uy*ux*cos_calc + uz*sin(radian);
            rotateMat[2] = uz*ux*cos_calc - uy*sin(radian);
            rotateMat[3] = 0;
            rotateMat[4] = ux*uy*cos_calc - uz*sin(radian);
            rotateMat[5] = cos(radian) + uy*uy*cos_calc;
            rotateMat[6] = uz*uy*cos_calc + ux*sin(radian);
            rotateMat[7] = 0;
            rotateMat[8] = ux*uz*cos_calc + uy*sin(radian);
            rotateMat[9] = uy*uz*cos_calc - ux*sin(radian);
            rotateMat[10] = cos(radian) + uz*uz*cos_calc;
            rotateMat[11] = 0;
            rotateMat[12] = 0;
            rotateMat[13] = 0;
            rotateMat[14] = 0;
            rotateMat[15] = 1;

            return rotateMat;
        }

        template<typename T>
        mat<4, 4, T> scale(T sx, T sy, T sz)
        {
            mat<4, 4, T> scaleMat;

            // TODO: Fill up this function properly 
            scaleMat[0] = sx;
            scaleMat[5] = sy;
            scaleMat[10] = sz;
            scaleMat[15] = 1;

            return scaleMat;
        }

        template<typename T>
        mat<4, 4, T> lookAt(T eyeX, T eyeY, T eyeZ, T centerX, T centerY, T centerZ, T upX, T upY, T upZ)
        {
            mat<4, 4, T> viewMat;

            // TODO: Fill up this function properly 
            T cam_pos[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, eyeX, eyeY, eyeZ, 1};
            T cam_x_axis, cam_y_axis, cam_z_axis, cam_x_axis_x, cam_x_axis_y, cam_x_axis_z, 
                cam_y_axis_x, cam_y_axis_y, cam_y_axis_z, cam_z_axis_x, cam_z_axis_y, cam_z_axis_z, 
                n_x, n_y, n_z, m_x, m_y, m_z;

            cam_z_axis = sqrt((centerX-eyeX)*(centerX-eyeX) + (centerY-eyeY)*(centerY-eyeY) + (centerZ-eyeZ)*(centerZ-eyeZ));
            cam_z_axis_x = (eyeX-centerX)/cam_z_axis;
            cam_z_axis_y = (eyeY-centerY)/cam_z_axis;
            cam_z_axis_z = (eyeZ-centerZ)/cam_z_axis;  //center - eye 는 안됨

            n_x = upY*cam_z_axis_z - upZ*cam_z_axis_y;
            n_y = upZ*cam_z_axis_x - upX*cam_z_axis_z;
            n_z = upX*cam_z_axis_y - upY*cam_z_axis_x;
            cam_x_axis = sqrt((n_x*n_x)+(n_y*n_y)+(n_z*n_z));
            cam_x_axis_x = n_x/cam_x_axis;
            cam_x_axis_y = n_y/cam_x_axis;
            cam_x_axis_z = n_z/cam_x_axis;

            m_x = cam_z_axis_y*cam_x_axis_z - cam_z_axis_z*cam_x_axis_y;
            m_y = cam_z_axis_z*cam_x_axis_x - cam_z_axis_x*cam_x_axis_z;
            m_z = cam_z_axis_x*cam_x_axis_y - cam_z_axis_y*cam_x_axis_x;
            cam_y_axis = sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
            cam_y_axis_x = m_x/cam_y_axis;
            cam_y_axis_y = m_y/cam_y_axis;
            cam_y_axis_z = m_z/cam_y_axis;

            T pos[16] = {cam_x_axis_x, cam_y_axis_x, cam_y_axis_x, 0, cam_x_axis_y, cam_y_axis_y, cam_y_axis_y, 0, cam_x_axis_z, cam_y_axis_z, cam_y_axis_z, 0, 0, 0, 0, 1};
            for(int i = 12; i < 15; i++){
                cam_pos[i] = (-1)*cam_pos[i];
            }

            int col = 0;
            while(col < 4){                
                int i = 0;
                for(int row = 0; row < 4; row++){
                    int n = 0;
                    int p = 0;
                    for(int m = 0; m < 4; m++){
                        viewMat[col+(4*i)] += pos[col+(4*n)]*cam_pos[4*row+p];
                        n++;
                        p++;
                    }
                    i++;       
                }
                col++;
            }
            

            return viewMat;
        }

        template<typename T>
        mat<4, 4, T> ortho(T left, T right, T bottom, T top, T nearVal, T farVal)
        {
            mat<4, 4, T> orthoMat;
            
            // TODO: Fill up this function properly 
            orthoMat[0] = 2/(right-left);
            orthoMat[5] = 2/(top-bottom);
            orthoMat[10] = (-1)*2/(farVal-nearVal);
            orthoMat[12] = (-1)*(right+1)/(right-1);
            orthoMat[13] = (-1)*(top+bottom)/(top-bottom);
            orthoMat[14] = (-1)*(farVal+nearVal)/(farVal-nearVal);
            orthoMat[15] = 1;

            return orthoMat;
        }

        template<typename T>
        mat<4, 4, T> frustum(T left, T right, T bottom, T top, T nearVal, T farVal)
        {
           mat<4, 4, T> frustumMat;

           // TODO: Fill up this function properly 
           frustumMat[0] = (2*nearVal)/(right-left);
           frustumMat[5] = (2*nearVal)/(top-bottom);
           frustumMat[8] = (right+left)/(right-left);
           frustumMat[9] = (top+bottom)/(top-bottom);
           frustumMat[10] = (-1)*(farVal+nearVal)/(farVal-nearVal);
           frustumMat[11] = -1;
           frustumMat[14] = (-1)*(2*farVal*nearVal)/(farVal-nearVal);

           return frustumMat;
        }

        template<typename T>
        mat<4, 4, T> perspective(T fovy, T aspect, T zNear, T zFar)
        {
          T  right = 0;
          T  top = 0;

          // TODO: Fill up this function properly 
          double radian = (M_PI/180)*fovy;
          
          top = zNear*tan(radian/2);
          right = aspect*top;

          return frustum(-right, right, -top, top, zNear, zFar);
        }
    }
}
#endif

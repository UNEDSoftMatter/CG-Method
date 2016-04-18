global_settings{ assumed_gamma 1.0 }
//#default{ finish{ ambient 0.1 diffuse 0.9 }} 
#default{ finish{ ambient 0.2 diffuse 0.7 }} 
//------------------------------------------------------------------------
#include "colors.inc"
#include "textures.inc"
#include "glass.inc"
#include "metals.inc"
#include "golds.inc"
#include "stones.inc"
#include "woods.inc"
#include "shapes.inc"
#include "shapes2.inc"
#include "functions.inc"
#include "math.inc"
#include "transforms.inc"
#include "objects.inc"

#declare Camera_1 = camera {/*ultra_wide_angle*/ angle 53 // diagonal view
                            //location  <0.7 , 0.6 ,-0.7>
                            // location  <0.5 , 0.6 ,-0.65>
                            //location  <20 , 30 ,30>
                            location  <60 , 90 ,90>
                            right     x*image_width/image_height
                            //look_at   < 0.15 , 0.0 , -0.135>}
                            look_at   < 0.15 , 0.15 , -0.135>}

camera{Camera_1}

// sun -------------------------------------------------------------------
light_source{<-1500,2500,2500> color White*0.9 }
// sky -------------------------------------------------------------------
sky_sphere{ pigment{ gradient <0,1,0>
                     color_map{ [0   color rgb<1,1,1>         ]//White
                                [0.4 color rgb<0.24,0.34,0.56>*0.8]//~Navy
                                [0.6 color rgb<0.24,0.34,0.56>*0.8]//~Navy
                                [1.0 color rgb<1,1,1>         ]//White
                              }
                   scale 2 }
           } // end of sky_sphere 


// ground -----------------------------------------------------------------
           //---------------------------------<<< settings of squared plane dimensions
           #declare RasterScale = 2.5;
           #declare RasterHalfLine  = 0.02;  
           #declare RasterHalfLineZ = 0.02; 
           //-------------------------------------------------------------------------
           #macro Raster(RScale, HLine) 
                  pigment{ gradient x scale RScale
                              color_map{    [0.000   color rgbt<1,1,1,0>*0.5]
                                            [0+HLine color rgbt<1,1,1,0>*0.5]
                                            [0+HLine color rgbt<1,1,1,1>]
                                            [1-HLine color rgbt<1,1,1,1>]
                                            [1-HLine color rgbt<1,1,1,0>*0.5]
                                            [1.000   color rgbt<1,1,1,0>*0.5]} }
         #end// of Raster(RScale, HLine)-macro    
//-------------------------------------------------------------------------
    
// squared plane XZ
plane { <0,1,0>, 0    // plane with layered textures
        texture { pigment{color SkyBlue*0.8}
                 normal {bumps 0.25 scale 0.75}
                 finish {ambient 0.1 diffuse 0.8}
                }
        texture { Raster(RasterScale,RasterHalfLine ) rotate<0,0,0> }
        texture { Raster(RasterScale,RasterHalfLineZ) rotate<0,90,0>}
        rotate<0,0,0>
      }
//------------------------------------------------ end of squared plane XZ

fog{ fog_type   2
     distance   10.0
     // color      rgb<1,1,1>*0.9
     color      SkyBlue*0.8 //rgb<0.1,0.1,0.9>
     fog_offset 1.0
     fog_alt    1.5
     turbulence 1.8
   } //---------------------------------------------------------------
 
//object {Boxes}
object {Particles}
object {Verlet}

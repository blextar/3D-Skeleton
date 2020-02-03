#include <cstdlib>
#include <list>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace cvlab
{

   struct point3d
   {
      double x, y, z;
      explicit point3d() : x(0), y(0), z(0){}
      point3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_){}
   };

   struct color_rgb24
   {
      unsigned char r, g, b;
      color_rgb24(unsigned char r_, unsigned char g_, unsigned char b_) : r(r_), g(g_), b(b_) {}
      explicit color_rgb24() : r(0), g(0), b(0) {}
   };

   const color_rgb24 white(255,255,255);
   const color_rgb24 black(0,0,0);
   const color_rgb24 light_grey(128,128,128);
   const color_rgb24 red(255,0,0);
   const color_rgb24 green(0,255,0);
   const color_rgb24 blue(0,0,255);
   const color_rgb24 magenta(255,0,255);
   const color_rgb24 cyan(0,255,255);
   const color_rgb24 yellow(255,255,0);
   const color_rgb24 pink(255,86,228);
   const color_rgb24 orange(207,88,0);
   const color_rgb24 forest(0,78,7);
   const color_rgb24 dark_grey(84,84,84);
   const color_rgb24 brown(139,69,19);
   const color_rgb24 mud(47,79,47);
   const color_rgb24 lime(124,252,0);
   const color_rgb24 ocra(139,101,8);
   const color_rgb24 cream(245,222,179);
   const color_rgb24 dark_red(128,0,0);

   class cloud3d
   {
      struct point3d_colored : public point3d
      {
         color_rgb24 color;
         explicit point3d_colored() : point3d() {}
         point3d_colored(const point3d& p, const color_rgb24& colr) : point3d(p), color(colr) {}
      };

      public:
         explicit cloud3d() : _npoints(0) {}

         void add_point(const point3d& p, const color_rgb24& color)
         {
            points.push_back(point3d_colored(p,color));
            ++_npoints;
         }

         void save_as_vtk(const std::string& fname) const;


      private:
         std::list<point3d_colored> points;
         unsigned int _npoints;
         unsigned int _nfaces;
   };

   void cloud3d::save_as_vtk(const std::string& fname) const
   {
      std::cout << "Saving " << fname << "..." << std::flush;

      std::ofstream model(fname.c_str(), std::ios::out);
      std::list<point3d_colored>::const_iterator it(points.begin()), it_end(points.end());

      model << "# vtk DataFile Version 2.0\ncvlab\nASCII\nDATASET UNSTRUCTURED_GRID" << std::endl;
      model << "POINTS " << _npoints << " float" << std::endl;

      // write xyz data
      while(it!=it_end) {
         model << std::fixed << std::setprecision(12) << it->x << " " << it->y << " " << it->z << std::endl;
         ++it;
      }

      model << "CELLS " << _npoints << " " << 2*_npoints << std::endl;

      // write cell data
      for (unsigned int i=0; i<_npoints; ++i)
         model << "1 " << i << "\n";

      // write cell types
      model << "CELL_TYPES " << _npoints << std::endl;
      for (unsigned int i=0; i<_npoints; ++i)
         model << "1 \n";

      // write z scalar values
      model << "\nPOINT_DATA " << _npoints << "\nVECTORS ScanColors unsigned_char" << std::endl;
      it = points.begin();
      while (it!=it_end){
         model << (unsigned int)it->color.r << " " << (unsigned int)it->color.g << " " << (unsigned int)it->color.b << "\n";
         ++it;
      }
      model << std::endl;

      model.close();
      std::cout << " done." << std::endl;
   }
}

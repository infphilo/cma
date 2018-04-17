/*
 * Copyright 2018, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of CMA.
 *
 * CMA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CMA.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include "tiffio.h"

const double PI = 3.141592653589793238463;

std::vector<std::vector<float> > net_mat_list;
std::vector<std::vector<float> > net_bias_list;
std::vector<std::pair<size_t, size_t> > net_dim_list;

bool init_network(const char* filename) {
  std::ifstream netfile(filename);
  if(!netfile.is_open())
    return false;

  size_t line = 0;
  netfile >> line;
  size_t cur_line = 0;
  while(netfile.good()) {
    std::vector<float> vec;
    size_t n = 1, m = 1;
    if(cur_line % 2 == 0) {
      netfile >> n >> m;
      net_dim_list.push_back(std::make_pair<size_t, size_t>(m, n));
    } else {
      netfile >> n;
    }    

    size_t mn = m * n;
    float number;
    while(true) {
      netfile >> number;
      vec.push_back(number);
      mn--;
      if(mn == 0) break;
    } 
    if(mn > 0) {
      std::cerr << "Error occurs in reading circle_mat.txt file" << std::endl;
    }

    if(cur_line % 2 == 0) {
      net_mat_list.push_back(vec);
    } else {
      net_bias_list.push_back(vec);
    }

    cur_line++;
    if(cur_line == line) {
      break;
    }
  }
  netfile.close();

  return true;
}

void draw_circle() {
  TIFF* tif = TIFFOpen("circle.tif", "w");
  if(tif == NULL) {
    return;
  }
  
  const uint32 width = 250, height = 250, radius = 100;
  const uint32 cx = width / 2, cy = height / 2;
  const uint32 nsamples = 4;
  
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, nsamples);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW) ;
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8) ;
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, width*nsamples));
  
  uint8* img = new uint8[height * width];
  memset(img, 0, width * height);
  
  const size_t num_theta = 1000;
  for(size_t i = 0; i < num_theta; i++) {
    float theta = 2 * PI * i / (num_theta - 1);
    float fx = cos(theta);
    float fy = sin(theta);
    uint32 x = cx + fx * radius;
    uint32 y = cy + fy * radius;
    assert(x >= 0 && x < width && y >= 0 && y < height);
    img[y * width + x] = 255;
  }
  
  for(uint32 r = 0; r < height; r++) {
    uint8* row = new uint8[width * nsamples];
    for(uint32 c = 0; c < width; c++) {
      row[c * nsamples] = img[r * width + c];
      row[c * nsamples + 1] = 0;
      row[c * nsamples + 2] = 0;
      row[c * nsamples + 3] = (row[c * nsamples] > 0 ? 255 : 0);
    }
    if(TIFFWriteScanline(tif, row, r) != 1) {
      std::cout<< "Unable to write a row." << std::endl;
      break;
    }
  }
  
  delete []img;
  TIFFClose(tif);
}

void get_xy(float theta,
	    float& x,
	    float& y) {
  assert(net_mat_list.size() == net_bias_list.size() && net_mat_list.size() == net_dim_list.size());
  std::vector<float> values;
  values.push_back(theta);
  std::vector<float> next_values;  
  for(size_t l = 0; l < net_dim_list.size(); l++) {
    size_t m = net_dim_list[l].first;
    size_t n = net_dim_list[l].second;

    // DK - debugging purposes
    // std::cout << "Layer: " << l << ", m: " << m << ", n: " << n << std::endl;
    
    for(size_t m2 = 0; m2 < m; m2++) {
      float value = 0.0f;
      for(size_t n2 = 0; n2 < n; n2++) {
	value += (net_mat_list[l][m2 + n2 * m] * values[n2]);
      }
      value += net_bias_list[l][m2];
      if(l + 1 < net_dim_list.size()) {
	if(value < 0) value = 0.0f;
      }

      // DK - debugging purposes
      // std::cout << "\tAt " << m2 << ", " << value << std::endl;
      
      next_values.push_back(value);
    }

    values = next_values;
    next_values.clear();
  }

  assert(values.size() == 2);
  x = values[0];
  y = values[1];
}

void draw_approx_circle() {
  if(!init_network("circle_mat_4L.txt")) {
    return;
  }
  
  TIFF* tif = TIFFOpen("circle_approx.tif", "w");
  if(tif == NULL) {
    return;
  }
  
  const uint32 width = 250, height = 250, radius = 100;
  const uint32 cx = width / 2, cy = height / 2;
  const uint32 nsamples = 1;
  
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW) ;
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8) ;
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, width*nsamples));
  
  uint8* img = new uint8[height * width];
  memset(img, 255, width * height);
  
  const size_t num_theta = 1000;
  for(size_t i = 0; i < num_theta; i++) {
    float theta = 2 * PI * i / (num_theta - 1);
    float fx = 0.0f, fy = 0.0f;
    get_xy(theta,
	   fx,
	   fy);
    uint32 x = cx + fx * radius;
    uint32 y = cy + fy * radius;
    assert(x >= 0 && x < width && y >= 0 && y < height);
    img[y * width + x] = 0;
  }

  // Draw tangent lines at 0, pi/4, pi/2, 2pi/3, pi, ...
#if 1
  for(size_t i = 0; i < 8; i++) {
    float theta = 2 * PI * i / 8;
    float fx = 0.0f, fy = 0.0f;
    get_xy(theta,
	   fx,
	   fy);
    uint32 x = cx + fx * radius;
    uint32 y = cy + fy * radius;
    float fx1 = 0.0f, fy1 = 0.0f;
    get_xy(theta - 0.01f,
	   fx1,
	   fy1);
    float fx2 = 0.0f, fy2 = 0.0f;
    get_xy(theta + 0.01f,
	   fx2,
	   fy2);
    float r = sqrt((fx2 - fx1) * (fx2 - fx1) + (fy2 - fy1) * (fy2 - fy1));
    float dx = (fx2 - fx1) / r;
    float dy = (fy2 - fy1) / r;
    for(int j = -50; j < 51; j++) {
      int32 ix = x + dx * j;
      if(ix < 0 || ix >= width) continue;
      int32 iy = y + dy * j;
      if(iy < 0 || iy >= height) continue;
      img[iy * width + ix] = 0;
    }
  }
#endif
  
  for(uint32 r = 0; r < height; r++) {
    if(TIFFWriteScanline(tif, &img[r * width], r) != 1) {
      std::cout<< "Unable to write a row." << std::endl;
      break;
    }
  }
  
  delete []img;
  TIFFClose(tif);
}

void draw_sphere() {
  TIFF* tif = TIFFOpen("sphere.tif", "w");
  if(tif == NULL) {
    return;
  }
  
  const uint32 width = 250, height = 250, depth = 250, radius = 100;
  const uint32 cx = width / 2, cy = height / 2, cz = depth / 2;
  const uint32 nsamples = 4;

  uint8* img = new uint8[depth * height * width];
  memset(img, 0, depth * width * height);
  
  const size_t num_theta = 1000, num_theta2 = 1000;
  for(size_t i = 0; i < num_theta; i++) {
    float theta = 2 * PI * i / (num_theta - 1);
    for(size_t j = 0; j < num_theta2; j++) {
      float theta2 = PI * j / (num_theta2 - 1);
      float fx = cos(theta) * sin(theta2);
      float fy = sin(theta) * sin(theta2);
      float fz = cos(theta2);
      uint32 x = cx + fx * radius;
      uint32 y = cy + fy * radius;
      uint32 z = cz + fz * radius;
      assert(x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth);
      img[(z * height + y) * width + x] = 255;
    }
  }

  // DK - debugging purposes
#if 0
  for(size_t z = 0; z < depth; z++) {
    if(z != 125) continue;
    for(size_t y = 0; y < height; y++) {
      for(size_t x = 0; x < width; x++) {
	uint8 v = img[(z * height + y) * width + x];
	std::cout << (v > 0 ? '1' : '0');
      }
      std::cout << std::endl;
    }

    break;
  }
#endif

  for(uint32 page = 0; page < depth; page++) {
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, nsamples);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, width * nsamples));
    
    TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
    TIFFSetField(tif, TIFFTAG_PAGENUMBER, page, depth);
  
    for(uint32 r = 0; r < height; r++) {
      uint8* row = new uint8[width * nsamples];
      for(uint32 c = 0; c < width; c++) {
	row[c * nsamples] = img[(page * height + r) * width + c];
	row[c * nsamples + 1] = 0;
	row[c * nsamples + 2] = 0;
	row[c * nsamples + 3] = (row[c * nsamples] > 0 ? 255 : 0);
      }
      if(TIFFWriteScanline(tif, row, r) != 1) {
	std::cout<< "Unable to write a row." << std::endl;
	break;
      }
    }
    TIFFWriteDirectory(tif);
  }  

  delete []img;
  TIFFClose(tif);
}

void get_xyz(float theta,
	     float theta2,
	     float& x,
	     float& y,
	     float& z) {
  assert(net_mat_list.size() == net_bias_list.size() && net_mat_list.size() == net_dim_list.size());
  std::vector<float> values;
  values.push_back(theta);
  values.push_back(theta2);
  std::vector<float> next_values;  
  for(size_t l = 0; l < net_dim_list.size(); l++) {
    size_t m = net_dim_list[l].first;
    size_t n = net_dim_list[l].second;

    // DK - debugging purposes
    // std::cout << "Layer: " << l << ", m: " << m << ", n: " << n << std::endl;
    
    for(size_t m2 = 0; m2 < m; m2++) {
      float value = 0.0f;
      for(size_t n2 = 0; n2 < n; n2++) {
	value += (net_mat_list[l][m2 + n2 * m] * values[n2]);
      }
      value += net_bias_list[l][m2];
      if(l + 1 < net_dim_list.size()) {
	if(value < 0) value = 0.0f;
      }

      // DK - debugging purposes
      // std::cout << "\tAt " << m2 << ", " << value << std::endl;
      
      next_values.push_back(value);
    }

    values = next_values;
    next_values.clear();
  }

  assert(values.size() == 3);
  x = values[0];
  y = values[1];
  z = values[2];
}

void draw_approx_sphere() {
  if(!init_network("sphere_mat_12L_3K.txt")) {
    return;
  }

  TIFF* tif = TIFFOpen("sphere_approx.tif", "w");
  if(tif == NULL) {
    return;
  }
  
  const uint32 width = 250, height = 250, depth = 250, radius = 100;
  const uint32 cx = width / 2, cy = height / 2, cz = depth / 2;
  const uint32 nsamples = 4;

  uint8* img = new uint8[depth * height * width];
  memset(img, 0, depth * width * height);
  
  const size_t num_theta = 1000, num_theta2 = 1000;
  for(size_t i = 0; i < num_theta; i++) {
    float theta = 2 * PI * i / (num_theta - 1);
    for(size_t j = 0; j < num_theta2; j++) {
      float theta2 = PI * j / (num_theta2 - 1);
      float fx = 0.0f, fy = 0.0f, fz = 0.0f;
      get_xyz(theta,
	      theta2,
	      fx,
	      fy,
	      fz);
      uint32 x = cx + fx * radius;
      uint32 y = cy + fy * radius;
      uint32 z = cz + fz * radius;
      assert(x >= 0 && x < width && y >= 0 && y < height && z >= 0 && z < depth);
      img[(z * height + y) * width + x] = 255;

      // DK - debugging purposes
      // std::cout << "x: " << x << ", y: " << y << ", z: " << z << std::endl;
    }

    // DK - debugging purposes
    if(i % 100 == 0) {
      std::cout << "i: " << i << std::endl;
    }
  }

    // DK - debugging purposes
#if 1
  for(size_t z = 0; z < depth; z++) {
    if(z != 125) continue;
    for(size_t y = 0; y < height; y++) {
      for(size_t x = 0; x < width; x++) {
	uint8 v = img[(z * height + y) * width + x];
	std::cout << (v > 0 ? '1' : '0');
      }
      std::cout << std::endl;
    }

    break;
  }
#endif


  // Draw tangent lines at 0, pi/4, pi/2, 2pi/3, pi, ...
#if 0
  for(size_t i = 0; i < 8; i++) {
    float theta = 2 * PI * i / 8;
    float fx = 0.0f, fy = 0.0f;
    get_xy(theta,
	   fx,
	   fy);
    uint32 x = cx + fx * radius;
    uint32 y = cy + fy * radius;
    float fx1 = 0.0f, fy1 = 0.0f;
    get_xy(theta - 0.01f,
	   fx1,
	   fy1);
    float fx2 = 0.0f, fy2 = 0.0f;
    get_xy(theta + 0.01f,
	   fx2,
	   fy2);
    float r = sqrt((fx2 - fx1) * (fx2 - fx1) + (fy2 - fy1) * (fy2 - fy1));
    float dx = (fx2 - fx1) / r;
    float dy = (fy2 - fy1) / r;
    for(int j = -50; j < 51; j++) {
      int32 ix = x + dx * j;
      if(ix < 0 || ix >= width) continue;
      int32 iy = y + dy * j;
      if(iy < 0 || iy >= height) continue;
      img[iy * width + ix] = 0;
    }
  }
#endif


  for(uint32 page = 0; page < depth; page++) {
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, nsamples);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, width * nsamples));
    
    TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
    TIFFSetField(tif, TIFFTAG_PAGENUMBER, page, depth);
  
    for(uint32 r = 0; r < height; r++) {
      uint8* row = new uint8[width * nsamples];
      for(uint32 c = 0; c < width; c++) {
	row[c * nsamples] = img[(page * height + r) * width + c];
	row[c * nsamples + 1] = 0;
	row[c * nsamples + 2] = 0;
	row[c * nsamples + 3] = (row[c * nsamples] > 0 ? 255 : 0);
      }
      if(TIFFWriteScanline(tif, row, r) != 1) {
	std::cout<< "Unable to write a row." << std::endl;
	break;
      }
    }
    TIFFWriteDirectory(tif);
  }  

  delete []img;
  TIFFClose(tif);
}

int main(int argc,char **argv) {
  // draw_circle();
  // draw_approx_circle();
  // draw_sphere();
  draw_approx_sphere();
  return 0;
}


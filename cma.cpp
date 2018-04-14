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

#include <stdio.h>
#include <stdlib.h>
#include "tiffio.h"

int main(int argc,char **argv)
{
   TIFF* tif = TIFFOpen("circle.tif", "w");
   if(tif == NULL) {
     return -1;
   }
   
   uint32 w, h;
   size_t npixels;
   uint32* raster;

   TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
   TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
   
   TIFFClose(tif);
   exit(0);
}


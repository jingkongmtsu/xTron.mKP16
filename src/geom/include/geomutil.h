#ifndef __GEOMUTIL_H__
#define __GEOMUTIL_H__
#include "libgen.h"

namespace geomutil {

	///
	/// counter-clockwize rotation matrix by X axis, 
	/// angle is in radian already
	///
	void rotationByX (Double& y, Double& z, Double angle);

	///
	/// counter-clockwize rotation matrix by Y axis, 
	/// angle is in radian already
	///
	void rotationByY (Double& x, Double& z, Double angle);

	///
	/// counter-clockwize rotation matrix by Z axis, 
	/// angle is in radian already
	///
	void rotationByZ (Double& x, Double& y, Double angle);

}

#endif

#ifndef GL_MATRIX_H
#define GL_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * gl-matrix.c - High performance matrix and vector operations for OpenGL
 * Version 1.2.3
 */
 
#define GL_MATRIX_MAJOR_VERSION 1
#define GL_MATRIX_MINOR_VERSION 2
#define GL_MATRIX_MICRO_VERSION 3

#define GL_MATRIX_VERSION  "1.2.3"

/* Hex version number. A value of 0x010203 means version 1.2.3.
  useful for comparisons. e.g. GL_MATRIX_VERSION_HEX >= 0x010203 */
#define GL_MATRIX_VERSION_HEX  ((GL_MATRIX_MAJOR_VERSION << 16) | \
                              (GL_MATRIX_MINOR_VERSION << 8) | \
                              (GL_MATRIX_MICRO_VERSION))

typedef double *ovr_vec3_t;
typedef double *ovr_vec4_t;
typedef double *ovr_mat3_t;
typedef double *ovr_mat4_t;
typedef double *quat_t;

/*
 * ovr_vec3_t - 3 Dimensional Vector
 */

/*
 * vec3_create
 * Creates a new instance of a ovr_vec3_t
 *
 * Params:
 * vec - Optional, ovr_vec3_t containing values to initialize with. If NULL, the 
 * vector will be initialized with zeroes.
 *
 * Returns:
 * New vec3
 */
ovr_vec3_t vec3_create(ovr_vec3_t vec);

/*
 * vec3_set
 * Copies the values of one ovr_vec3_t to another
 *
 * Params:
 * vec - ovr_vec3_t containing values to copy
 * dest - ovr_vec3_t receiving copied values
 *
 * Returns:
 * dest
 */
ovr_vec3_t vec3_set(ovr_vec3_t vec, ovr_vec3_t dest);

/*
 * vec3_add
 * Performs a vector addition
 *
 * Params:
 * vec - vec3, first operand
 * vec2 - vec3, second operand
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_add(ovr_vec3_t vec, ovr_vec3_t vec2, ovr_vec3_t dest);

/*
 * vec3_subtract
 * Performs a vector subtraction
 *
 * Params:
 * vec - vec3, first operand
 * vec2 - vec3, second operand
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_subtract(ovr_vec3_t vec, ovr_vec3_t vec2, ovr_vec3_t dest);

/*
 * vec3_multiply
 * Performs a vector multiplication
 *
 * Params:
 * vec - vec3, first operand
 * vec2 - vec3, second operand
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_multiply(ovr_vec3_t vec, ovr_vec3_t vec2, ovr_vec3_t dest);

/*
 * vec3_negate
 * Negates the components of a vec3
 *
 * Params:
 * vec - ovr_vec3_t to negate
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_negate(ovr_vec3_t vec, ovr_vec3_t dest);

/*
 * vec3_scale
 * Multiplies the components of a ovr_vec3_t by a scalar value
 *
 * Params:
 * vec - ovr_vec3_t to scale
 * val - Numeric value to scale by
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_scale(ovr_vec3_t vec, double val, ovr_vec3_t dest);

/*
 * vec3_normalize
 * Generates a unit vector of the same direction as the provided vec3
 * If vector length is 0, returns [0, 0, 0]
 *
 * Params:
 * vec - ovr_vec3_t to normalize
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_normalize(ovr_vec3_t vec, ovr_vec3_t dest);

/*
 * vec3_cross
 * Generates the cross product of two vec3s
 *
 * Params:
 * vec - vec3, first operand
 * vec2 - vec3, second operand
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_cross (ovr_vec3_t vec, ovr_vec3_t vec2, ovr_vec3_t dest);

/*
 * vec3_length
 * Caclulates the length of a vec3
 *
 * Params:
 * vec - ovr_vec3_t to calculate length of
 *
 * Returns:
 * Length of vec
 */
double vec3_length(ovr_vec3_t vec);

/*
 * vec3_dot
 * Caclulates the dot product of two vec3s
 *
 * Params:
 * vec - vec3, first operand
 * vec2 - vec3, second operand
 *
 * Returns:
 * Dot product of vec and vec2
 */
double vec3_dot(ovr_vec3_t vec, ovr_vec3_t vec2);

/*
 * vec3_direction
 * Generates a unit vector pointing from one vector to another
 *
 * Params:
 * vec - origin vec3
 * vec2 - ovr_vec3_t to point to
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_direction (ovr_vec3_t vec, ovr_vec3_t vec2, ovr_vec3_t dest);

/*
 * vec3_lerp
 * Performs a linear interpolation between two vec3
 *
 * Params:
 * vec - vec3, first vector
 * vec2 - vec3, second vector
 * lerp - interpolation amount between the two inputs
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */

ovr_vec3_t vec3_lerp(ovr_vec3_t vec, ovr_vec3_t vec2, double lerp, ovr_vec3_t dest);

/*
 * vec3_dist
 * Calculates the euclidian distance between two vec3
 *
 * Params:
 * vec - vec3, first vector
 * vec2 - vec3, second vector
 *
 * Returns:
 * distance between vec and vec2
 */
double vec3_dist(ovr_vec3_t vec, ovr_vec3_t vec2);

/*
 * vec3_unproject
 * Projects the specified ovr_vec3_t from screen space into object space
 * Based on Mesa gluUnProject implementation at: 
 * http://webcvs.freedesktop.org/mesa/Mesa/src/glu/mesa/project.c?revision=1.4&view=markup
 *
 * Params:
 * vec - vec3, screen-space vector to project
 * view - mat4, View matrix
 * proj - mat4, Projection matrix
 * viewport - vec4, Viewport as given to gl.viewport [x, y, width, height]
 * dest - Optional, ovr_vec3_t receiving unprojected result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t vec3_unproject(ovr_vec3_t vec, ovr_mat4_t view, ovr_mat4_t proj, ovr_vec4_t viewport, ovr_vec3_t dest);

/*
 * vec3_str
 * Writes a string representation of a vector
 *
 * Params:
 * vec - ovr_vec3_t to represent as a string
 * buffer - char * to store the results
 */
void vec3_str(ovr_vec3_t vec, char *buffer);

/*
 * ovr_mat3_t - 3x3 Matrix
 */

/*
 * mat3_create
 * Creates a new instance of a ovr_mat3_t
 *
 * Params:
 * mat - Optional, ovr_mat3_t containing values to initialize with. If NULL the result
 * will be initialized with zeroes.
 *
 * Returns:
 * New mat3
 */
ovr_mat3_t mat3_create(ovr_mat3_t mat);

/*
 * mat3_set
 * Copies the values of one ovr_mat3_t to another
 *
 * Params:
 * mat - ovr_mat3_t containing values to copy
 * dest - ovr_mat3_t receiving copied values
 *
 * Returns:
 * dest
 */
ovr_mat3_t mat3_set(ovr_mat3_t mat, ovr_mat3_t dest);

/*
 * mat3_identity
 * Sets a ovr_mat3_t to an identity matrix
 *
 * Params:
 * dest - ovr_mat3_t to set
 *
 * Returns:
 * dest
 */
 ovr_mat3_t mat3_identity(ovr_mat3_t dest);

/*
 * mat4.transpose
 * Transposes a ovr_mat3_t (flips the values over the diagonal)
 *
 * Params:
 * mat - ovr_mat3_t to transpose
 * dest - Optional, ovr_mat3_t receiving transposed values. If NULL, result is written to mat
 *
 * Returns:
 * dest is specified, mat otherwise
 */
ovr_mat3_t ovr_mat3_transpose(ovr_mat3_t mat, ovr_mat3_t dest);

/*
 * ovr_mat3_toMat4
 * Copies the elements of a ovr_mat3_t into the upper 3x3 elements of a mat4
 *
 * Params:
 * mat - ovr_mat3_t containing values to copy
 * dest - Optional, ovr_mat4_t receiving copied values
 *
 * Returns:
 * dest if not NULL, a new ovr_mat4_t otherwise
 */
ovr_mat4_t ovr_mat3_toMat4(ovr_mat3_t mat, ovr_mat4_t dest);

/*
 * mat3_str
 * Writes a string representation of a mat3
 *
 * Params:
 * mat - ovr_mat3_t to represent as a string
 * buffer - char * to store the results 
 */
void mat3_str(ovr_mat3_t mat, char *buffer);

/*
 * ovr_mat4_t - 4x4 Matrix
 */

/*
 * mat4_create
 * Creates a new instance of a ovr_mat4_t
 *
 * Params:
 * mat - Optional, ovr_mat4_t containing values to initialize with
 *
 * Returns:
 * New mat4
 */
ovr_mat4_t mat4_create(ovr_mat4_t mat);

/*
 * mat4_set
 * Copies the values of one ovr_mat4_t to another
 *
 * Params:
 * mat - ovr_mat4_t containing values to copy
 * dest - ovr_mat4_t receiving copied values
 *
 * Returns:
 * dest
 */
ovr_mat4_t mat4_set(ovr_mat4_t mat, ovr_mat4_t dest);

/*
 * mat4_identity
 * Sets a ovr_mat4_t to an identity matrix
 *
 * Params:
 * dest - ovr_mat4_t to set
 *
 * Returns:
 * dest
 */
ovr_mat4_t mat4_identity(ovr_mat4_t dest);
     
/*
 * ovr_mat4_transpose
 * Transposes a ovr_mat4_t (flips the values over the diagonal)
 *
 * Params:
 * mat - ovr_mat4_t to transpose
 * dest - Optional, ovr_mat4_t receiving transposed values. If NULL, result is written to mat
 *
 * Returns:
 * dest is specified, mat otherwise
 */
ovr_mat4_t ovr_mat4_transpose(ovr_mat4_t mat, ovr_mat4_t dest);

/*
 * mat4_determinant
 * Calculates the determinant of a mat4
 *
 * Params:
 * mat - ovr_mat4_t to calculate determinant of
 *
 * Returns:
 * determinant of mat
 */
double mat4_determinant(ovr_mat4_t mat);

/*
 * mat4_inverse
 * Calculates the inverse matrix of a mat4
 *
 * Params:
 * mat - ovr_mat4_t to calculate inverse of
 * dest - Optional, ovr_mat4_t receiving inverse matrix. If NULL, result is written to mat
 *
 * Returns:
 * dest is specified, mat otherwise, NULL if matrix cannot be inverted
 */
ovr_mat4_t mat4_inverse(ovr_mat4_t mat, ovr_mat4_t dest);

/*
 * ovr_mat4_toRotationMat
 * Copies the upper 3x3 elements of a ovr_mat4_t into another mat4
 *
 * Params:
 * mat - ovr_mat4_t containing values to copy
 * dest - Optional, ovr_mat4_t receiving copied values
 *
 * Returns:
 * dest is specified, a new ovr_mat4_t otherwise
 */
ovr_mat4_t ovr_mat4_toRotationMat(ovr_mat4_t mat, ovr_mat4_t dest);

/*
 * ovr_mat4_toMat3
 * Copies the upper 3x3 elements of a ovr_mat4_t into a mat3
 *
 * Params:
 * mat - ovr_mat4_t containing values to copy
 * dest - Optional, ovr_mat3_t receiving copied values
 *
 * Returns:
 * dest is specified, a new ovr_mat3_t otherwise
 */
ovr_mat3_t ovr_mat4_toMat3(ovr_mat4_t mat, ovr_mat3_t dest);

/*
 * ovr_mat4_toInverseMat3
 * Calculates the inverse of the upper 3x3 elements of a ovr_mat4_t and copies the result into a mat3
 * The resulting matrix is useful for calculating transformed normals
 *
 * Params:
 * mat - ovr_mat4_t containing values to invert and copy
 * dest - Optional, ovr_mat3_t receiving values
 *
 * Returns:
 * dest is specified, a new ovr_mat3_t otherwise, NULL if the matrix cannot be inverted
 */
ovr_mat3_t ovr_mat4_toInverseMat3(ovr_mat4_t mat, ovr_mat3_t dest);

/*
 * mat4_multiply
 * Performs a matrix multiplication
 *
 * Params:
 * mat - mat4, first operand
 * mat2 - mat4, second operand
 * dest - Optional, ovr_mat4_t receiving operation result. If NULL, result is written to mat
 *
 * Returns:
 * dest if not NULL, mat otherwise
 */
ovr_mat4_t mat4_multiply(ovr_mat4_t mat, ovr_mat4_t mat2, ovr_mat4_t dest);

/*
 * mat4_multiplyVec3
 * Transforms a ovr_vec3_t with the given matrix
 * 4th vector component is implicitly '1'
 *
 * Params:
 * mat - ovr_mat4_t to transform the vector with
 * vec - ovr_vec3_t to transform
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_mat4_t mat4_multiplyVec3(ovr_mat4_t mat, ovr_vec3_t vec, ovr_mat4_t dest);

/*
 * mat4_multiplyVec4
 * Transforms a vec4 with the given matrix
 *
 * Params:
 * mat - ovr_mat4_t to transform the vector with
 * vec - vec4 to transform
 * dest - Optional, vec4 receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_mat4_t mat4_multiplyVec4(ovr_mat4_t mat, ovr_vec4_t vec, ovr_mat4_t dest);

/*
 * ovr_mat4_translate
 * Translates a matrix by the given vector
 *
 * Params:
 * mat - ovr_mat4_t to translate
 * vec - ovr_vec3_t specifying the translation
 * dest - Optional, ovr_mat4_t receiving operation result. If NULL, result is written to mat
 *
 * Returns:
 * dest if not NULL, mat otherwise
 */
ovr_mat4_t ovr_mat4_translate(ovr_mat4_t mat, ovr_vec3_t vec, ovr_mat4_t dest);

/*
 * mat4_scale
 * Scales a matrix by the given vector
 *
 * Params:
 * mat - ovr_mat4_t to scale
 * vec - ovr_vec3_t specifying the scale for each axis
 * dest - Optional, ovr_mat4_t receiving operation result. If NULL, result is written to mat
 *
 * Returns:
 * dest if not NULL, mat otherwise
 */
ovr_mat4_t mat4_scale(ovr_mat4_t mat, ovr_vec3_t vec, ovr_mat4_t dest);

/*
 * mat4_rotate
 * Rotates a matrix by the given angle around the specified axis
 * If rotating around a primary axis (X,Y,Z) one of the specialized rotation functions should be used instead for performance
 *
 * Params:
 * mat - ovr_mat4_t to rotate
 * angle - angle (in radians) to rotate
 * axis - ovr_vec3_t representing the axis to rotate around 
 * dest - Optional, ovr_mat4_t receiving operation result. If NULL, result is written to mat
 *
 * Returns:
 * dest if not NULL, mat otherwise
 */
ovr_mat4_t mat4_rotate(ovr_mat4_t mat, double angle, ovr_vec3_t axis, ovr_mat4_t dest);

/*
 * mat4_rotateX
 * Rotates a matrix by the given angle around the X axis
 *
 * Params:
 * mat - ovr_mat4_t to rotate
 * angle - angle (in radians) to rotate
 * dest - Optional, ovr_mat4_t receiving operation result. If NULL, result is written to mat
 *
 * Returns:
 * dest if not NULL, mat otherwise
 */
ovr_mat4_t mat4_rotateX(ovr_mat4_t mat, double angle, ovr_mat4_t dest);

/*
 * mat4_rotateY
 * Rotates a matrix by the given angle around the Y axis
 *
 * Params:
 * mat - ovr_mat4_t to rotate
 * angle - angle (in radians) to rotate
 * dest - Optional, ovr_mat4_t receiving operation result. If NULL, result is written to mat
 *
 * Returns:
 * dest if not NULL, mat otherwise
 */
ovr_mat4_t mat4_rotateY(ovr_mat4_t mat, double angle, ovr_mat4_t dest);

/*
 * mat4_rotateZ
 * Rotates a matrix by the given angle around the Z axis
 *
 * Params:
 * mat - ovr_mat4_t to rotate
 * angle - angle (in radians) to rotate
 * dest - Optional, ovr_mat4_t receiving operation result. If NULL, result is written to mat
 *
 * Returns:
 * dest if not NULL, mat otherwise
 */
ovr_mat4_t mat4_rotateZ(ovr_mat4_t mat, double angle, ovr_mat4_t dest);

/*
 * mat4_frustum
 * Generates a frustum matrix with the given bounds
 *
 * Params:
 * left, right - scalar, left and right bounds of the frustum
 * bottom, top - scalar, bottom and top bounds of the frustum
 * near, far - scalar, near and far bounds of the frustum
 * dest - Optional, ovr_mat4_t frustum matrix will be written into
 *
 * Returns:
 * dest if not NULL, a new ovr_mat4_t otherwise
 */
ovr_mat4_t mat4_frustum(double left, double right, double bottom, double top, double near, double far, ovr_mat4_t dest);

/*
 * mat4_perspective
 * Generates a perspective projection matrix with the given bounds
 *
 * Params:
 * fovy - scalar, vertical field of view
 * aspect - scalar, aspect ratio. typically viewport width/height
 * near, far - scalar, near and far bounds of the frustum
 * dest - Optional, ovr_mat4_t frustum matrix will be written into
 *
 * Returns:
 * dest if not NULL, a new ovr_mat4_t otherwise
 */
ovr_mat4_t mat4_perspective(double fovy, double aspect, double near, double far, ovr_mat4_t dest);

/*
 * mat4_ortho
 * Generates a orthogonal projection matrix with the given bounds
 *
 * Params:
 * left, right - scalar, left and right bounds of the frustum
 * bottom, top - scalar, bottom and top bounds of the frustum
 * near, far - scalar, near and far bounds of the frustum
 * dest - Optional, ovr_mat4_t frustum matrix will be written into
 *
 * Returns:
 * dest if not NULL, a new ovr_mat4_t otherwise
 */
ovr_mat4_t mat4_ortho(double left, double right, double bottom, double top, double near, double far, ovr_mat4_t dest);

/*
 * mat4_lookAt
 * Generates a look-at matrix with the given eye position, focal point, and up axis
 *
 * Params:
 * eye - vec3, position of the viewer
 * center - vec3, point the viewer is looking at
 * up - ovr_vec3_t pointing "up"
 * dest - Optional, ovr_mat4_t frustum matrix will be written into
 *
 * Returns:
 * dest if not NULL, a new ovr_mat4_t otherwise
 */
ovr_mat4_t mat4_lookAt(ovr_vec3_t eye, ovr_vec3_t center, ovr_vec3_t up, ovr_mat4_t dest);

/*
 * mat4_fromRotationTranslation
 * Creates a matrix from a quaternion rotation and vector translation
 * This is equivalent to (but much faster than):
 *
 *     mat4_identity(dest);
 *     ovr_mat4_translate(dest, vec);
 *     ovr_mat4_t quatMat = mat4_create();
 *     quat_toMat4(quat, quatMat);
 *     mat4_multiply(dest, quatMat);
 *
 * Params:
 * quat - quat specifying the rotation by
 * vec - ovr_vec3_t specifying the translation
 * dest - Optional, ovr_mat4_t receiving operation result. If NULL, result is written to a new mat4
 *
 * Returns:
 * dest if not NULL, a new ovr_mat4_t otherwise
 */
ovr_mat4_t mat4_fromRotationTranslation(quat_t quat, ovr_vec3_t vec, ovr_mat4_t dest);

/*
 * mat4_str
 * Writes a string representation of a mat4
 *
 * Params:
 * mat - ovr_mat4_t to represent as a string
 * buffer - char * to store the results
 */
void mat4_str(ovr_mat4_t mat, char *buffer);

/*
 * quat - Quaternions 
 */

/*
 * quat_create
 * Creates a new instance of a quat_t
 *
 * Params:
 * quat - Optional, quat_t containing values to initialize with
 *
 * Returns:
 * New quat_t
 */
quat_t quat_create(quat_t quat);

/*
 * quat_set
 * Copies the values of one quat_t to another
 *
 * Params:
 * quat - quat_t containing values to copy
 * dest - quat_t receiving copied values
 *
 * Returns:
 * dest
 */
quat_t quat_set(quat_t quat, quat_t dest);

/*
 * quat_calculateW
 * Calculates the W component of a quat_t from the X, Y, and Z components.
 * Assumes that quaternion is 1 unit in length. 
 * Any existing W component will be ignored. 
 *
 * Params:
 * quat - quat_t to calculate W component of
 * dest - Optional, quat_t receiving calculated values. If NULL, result is written to quat
 *
 * Returns:
 * dest if not NULL, quat otherwise
 */
quat_t quat_calculateW(quat_t quat, quat_t dest);

/**
 * quat_dot
 * Calculates the dot product of two quaternions
 *
 * @param {quat4} quat First operand
 * @param {quat4} quat2 Second operand
 *
 * @return {number} Dot product of quat and quat2
 */
double quat_dot(quat_t quat, quat_t quat2);

/*
 * quat_inverse
 * Calculates the inverse of a quat_t
 *
 * Params:
 * quat - quat_t to calculate inverse of
 * dest - Optional, quat_t receiving inverse values. If NULL, result is written to quat
 *
 * Returns:
 * dest if not NULL, quat otherwise
 */
quat_t quat_inverse(quat_t quat, quat_t dest);

/*
 * quat_conjugate
 * Calculates the conjugate of a quat_t
 *
 * Params:
 * quat - quat_t to calculate conjugate of
 * dest - Optional, quat_t receiving conjugate values. If NULL, result is written to quat
 *
 * Returns:
 * dest if not NULL, quat otherwise
 */
quat_t quat_conjugate(quat_t quat, quat_t dest);

/*
 * quat_length
 * Calculates the length of a quat_t
 *
 * Params:
 * quat - quat_t to calculate length of
 *
 * Returns:
 * Length of quat
 */
double quat_length(quat_t quat);

/*
 * quat_normalize
 * Generates a unit quaternion of the same direction as the provided quat_t
 * If quaternion length is 0, returns [0, 0, 0, 0]
 *
 * Params:
 * quat - quat_t to normalize
 * dest - Optional, quat_ receiving operation result. If NULL, result is written to quat
 *
 * Returns:
 * dest if not NULL, quat otherwise
 */
quat_t quat_normalize(quat_t quat, quat_t dest);

/*
 * quat_multiply
 * Performs a quaternion multiplication
 *
 * Params:
 * quat - quat_t, first operand
 * quat2 - quat_t, second operand
 * dest - Optional, quat_t receiving operation result. If NULL, result is written to quat
 *
 * Returns:
 * dest if not NULL, quat otherwise
 */
quat_t quat_multiply(quat_t quat, quat_t quat2, quat_t dest);

/*
 * quat_multiplyVec3
 * Transforms a ovr_vec3_t with the given quaternion
 *
 * Params:
 * quat - quat_t to transform the vector with
 * vec - ovr_vec3_t to transform
 * dest - Optional, ovr_vec3_t receiving operation result. If NULL, result is written to vec
 *
 * Returns:
 * dest if not NULL, vec otherwise
 */
ovr_vec3_t quat_multiplyVec3(quat_t quat, ovr_vec3_t vec, ovr_vec3_t dest);

/*
 * quat_toMat3
 * Calculates a 3x3 matrix from the given quat_t
 *
 * Params:
 * quat - quat_t to create matrix from
 * dest - Optional, ovr_mat3_t receiving operation result
 *
 * Returns:
 * dest if not NULL, a new ovr_mat3_t otherwise
 */
ovr_mat3_t quat_toMat3(quat_t quat, ovr_mat3_t dest);

/*
 * quat_toMat4
 * Calculates a 4x4 matrix from the given quat_t
 *
 * Params:
 * quat - quat_t to create matrix from
 * dest - Optional, ovr_mat4_t receiving operation result
 *
 * Returns:
 * dest if not NULL, a new ovr_mat4_t otherwise
 */
quat_t quat_toMat4(quat_t quat, ovr_mat4_t dest);

/*
 * quat_slerp
 * Performs a spherical linear interpolation between two quat_t
 *
 * Params:
 * quat - quat_t, first quaternion
 * quat2 - quat_t, second quaternion
 * slerp - interpolation amount between the two inputs
 * dest - Optional, quat_t receiving operation result. If NULL, result is written to quat
 *
 * Returns:
 * dest if not NULL, quat otherwise
 */
quat_t quat_slerp(quat_t quat, quat_t quat2, double slerp, quat_t dest);

/*
 * quat_str
 * Writes a string representation of a quaternion
 *
 * Params:
 * quat - quat_t to represent as a string
 * buffer - char * to store the results
 */
void quat_str(quat_t quat, char *buffer);

#ifdef __cplusplus
}
#endif

#endif

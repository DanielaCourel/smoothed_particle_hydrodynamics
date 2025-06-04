#version 330 core // Use OpenGL 3.3 core profile

layout (location = 0) in vec3 aPos; // Input vertex position from VBO

uniform mat4 u_projection;  // Projection matrix from host
uniform mat4 u_modelview;   // Modelview matrix from host
uniform float u_particleSize; // Particle size from host (via ImGui)

void main() {
    // Transform position by modelview and projection matrices
    gl_Position = u_projection * u_modelview * vec4(aPos, 1.0);

    // Set the point size for gl_Points primitive.
    // gl_PointSize is a built-in output variable in vertex shader.
    gl_PointSize = u_particleSize;
}
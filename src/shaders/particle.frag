#version 330 core

out vec4 FragColor; // Output fragment color

void main() {
    // gl_PointCoord is a built-in input variable for point primitives.
    // It ranges from (0,0) to (1,1) across the point.
    // Calculate distance from center (0.5, 0.5) to create a circular point.
    float dist = distance(gl_PointCoord, vec2(0.5, 0.5));

    if (dist > 0.5) {
        // Discard fragments outside the circle to make points look round
        discard;
    } else {
        // Simple white color for particles
        FragColor = vec4(1.0, 1.0, 1.0, 1.0);
    }
}
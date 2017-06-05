#version 430 core

//in vec3 vertexPosition;
layout (location = 0) in vec3 aPos;

void main() {
        gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0f);
}

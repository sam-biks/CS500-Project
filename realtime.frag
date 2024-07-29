#version 330
uniform int mode;
in vec3 normalVec, lightVec, worldPos;
in vec2 texCoord;
out vec4 FragColor;

uniform bool emitter;
uniform vec3 ambient;
uniform vec3 Kd;
uniform vec3 Ks;
uniform float alpha;
uniform sampler2D tex;
uniform int lightNum;
uniform vec3 lightPosn[8];
uniform vec3 lightEmit[8];
uniform vec3 eyePos;
void main()
{
    if (emitter) {
        FragColor.xyz = Kd;
        return; }

    vec3 tKd = Kd;
    if (textureSize(tex,0).x>1)
        tKd = texture(tex, texCoord).xyz;

    vec3 N = normalize(normalVec);
    vec3 V = normalize(eyePos - worldPos);
    FragColor.xyz = ambient*tKd;
    for (int i=0;  i<lightNum;  i++) {
      vec3 L = normalize(lightPosn[i]-worldPos);
      vec3 H = normalize(L+V);

      float NL = max(0.0, dot(N,L));
      float HN = max(0.0, dot(H,N));

      FragColor.xyz += lightEmit[i]*(tKd*NL/3.14 + alpha/6.28*Ks*pow(HN,alpha)); }
};

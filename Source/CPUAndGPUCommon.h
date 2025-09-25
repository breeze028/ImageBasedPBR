#pragma once

#ifdef __cplusplus
#include "DirectXMath/DirectXMath.h"
typedef XMFLOAT4X4 float4x4;
typedef XMFLOAT4X3 float4x3;
typedef XMFLOAT2 float2;
typedef XMFLOAT3 float3;
typedef XMFLOAT4 float4;
#endif

#ifdef __cplusplus
#define SALIGN alignas(256)
#else
#define SALIGN
#endif

#define IBL_MODE_SPLIT_SUM_APPROXIMATION 0
#define IBL_MODE_SPHERICAL_HARMONICS_EXPONENTIAL 1
#define IBL_MODE_REFERENCE 2

#define MATERIAL_MODE_DIFFUSE_AND_SPECULAR 0
#define MATERIAL_MODE_DIFFUSE_ONLY 1
#define MATERIAL_MODE_SPECULAR_ONLY 2

struct SALIGN FPerDrawConstantData
{
	float4x4 ObjectToClip;
	float4x3 ObjectToWorld;
	float3 Albedo;
	float Metallic;
	float Roughness;
	float AO;
};

struct SALIGN FPerFrameConstantData
{
	float4 LightPositions[4];
	float4 LightColors[4];
	float4 ViewerPosition;
	float2 ViewportSize;
	int MaterialMode;
	int IBLMode;
	uint32_t NumFrames;
	bool bUseProgressiveRendering;
};

#ifdef __cplusplus
#undef SALIGN
#endif

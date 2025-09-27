#include "../CPUAndGPUCommon.h"
#include "Common.hlsli"

#define GRootSignature                                                                                     \
	"RootFlags(ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT), "                                                      \
	"CBV(b0), "                                                                                            \
	"DescriptorTable(CBV(b1), SRV(t0), SRV(t1), SRV(t2), SRV(t3), SRV(t4), visibility = SHADER_VISIBILITY_PIXEL), " \
	"StaticSampler("                                                                                       \
	"s0, "                                                                                                 \
	"filter = FILTER_MIN_MAG_MIP_LINEAR, "                                                                 \
	"visibility = SHADER_VISIBILITY_PIXEL, "                                                               \
	"addressU = TEXTURE_ADDRESS_BORDER, "                                                                  \
	"addressV = TEXTURE_ADDRESS_BORDER, "                                                                  \
	"addressW = TEXTURE_ADDRESS_BORDER)"

ConstantBuffer<FPerDrawConstantData> GPerDrawCB : register(b0);
ConstantBuffer<FPerFrameConstantData> GPerFrameCB : register(b1);
TextureCube GIrradianceMap : register(t0);
TextureCube GPrefilteredEnvMap : register(t1);
TextureCube GEnvMap : register(t2);
Texture2D GBRDFIntegrationMap : register(t3);
Texture2D GAccumulationBuffer : register(t4);
SamplerState GSampler : register(s0);

// Trowbridge-Reitz GGX normal distribution function.
float DistributionGGX(float3 N, float3 H, float Roughness)
{
	float Alpha = Roughness * Roughness;
	float Alpha2 = Alpha * Alpha;
	float NoH = dot(N, H);
	float NoH2 = NoH * NoH;
	float K = NoH2 * Alpha2 + (1.0f - NoH2);
	return Alpha2 / (PI * K * K);
}

float3 FresnelSchlick(float CosTheta, float3 F0)
{
	return F0 + (1.0f - F0) * pow(1.0f - CosTheta, 5.0f);
}

float3 FresnelSchlickRoughness(float CosTheta, float3 F0, float Roughness)
{
	return F0 + (max(1.0f - Roughness, F0) - F0) * pow(1.0f - CosTheta, 5.0f);
}

[RootSignature(GRootSignature)] void MainVS(
	in float3 InPosition : _Position,
	in float3 InNormal : _Normal,
	out float4 OutPosition : SV_Position,
	out float3 OutPositionWS : _Position,
	out float3 OutNormalWS : _Normal)
{
	OutPosition = mul(float4(InPosition, 1.0f), GPerDrawCB.ObjectToClip);
	OutPositionWS = mul(float4(InPosition, 1.0f), GPerDrawCB.ObjectToWorld);
	OutNormalWS = mul(InNormal, (float3x3)GPerDrawCB.ObjectToWorld);
}

	[RootSignature(GRootSignature)] void MainPS(
		in float4 InPosition : SV_Position,
		in float3 InPositionWS : _Position,
		in float3 InNormalWS : _Normal,
		out float4 OutColor : SV_Target0)
{
	float3 V = normalize(GPerFrameCB.ViewerPosition.xyz - InPositionWS);
	float3 N = normalize(InNormalWS);
	float NoV = saturate(dot(N, V));

	float3 Albedo = GPerDrawCB.Albedo;
	float Roughness = GPerDrawCB.Roughness;
	float Metallic = GPerDrawCB.Metallic;
	float AO = GPerDrawCB.AO;

	float3 F0 = float3(0.04f, 0.04f, 0.04f);
	F0 = lerp(F0, Albedo, Metallic);

	float3 Lo = 0.0f;
	// We manually disable point lights, only use image based lighting.
	for (int LightIdx = 0; LightIdx < 0; ++LightIdx)
	{
		float3 LightVector = GPerFrameCB.LightPositions[LightIdx].xyz - InPositionWS;

		float3 L = normalize(LightVector);
		float3 H = normalize(L + V);
		float NoL = saturate(dot(N, L));
		float HoV = saturate(dot(H, V));

		float Attenuation = 1.0f / dot(LightVector, LightVector);
		float3 Radiance = GPerFrameCB.LightColors[LightIdx].rgb * Attenuation;

		float3 F = FresnelSchlick(HoV, F0);

		float NDF = DistributionGGX(N, H, Roughness);
		float G = GeometrySmith(NoL, NoV, (Roughness + 1.0f) * 0.5f);

		float3 Specular = (NDF * G * F) / max(4.0f * NoV * NoL, 0.001f);

		float3 KS = F;
		float3 KD = 1.0f - KS;
		KD *= 1.0f - Metallic;

		Lo += (KD * (Albedo / PI) + Specular) * Radiance * NoL;
	}

	float3 R = reflect(-V, N);

	float3 F = FresnelSchlickRoughness(NoV, F0, Roughness);

	float3 KD = 1.0f - F;
	KD *= 1.0f - Metallic;

	float3 Diffuse = 0.0f;
	float3 Specular = 0.0f;
	if (GPerFrameCB.IBLMode == IBL_MODE_SPLIT_SUM_APPROXIMATION)
	{
		float3 Irradiance = GIrradianceMap.SampleLevel(GSampler, N, 0.0f).rgb;
		Diffuse = Irradiance * Albedo;

		float3 PrefilteredColor = GPrefilteredEnvMap.SampleLevel(GSampler, R, Roughness * 5.0f).rgb;
		float2 EnvBRDF = GBRDFIntegrationMap.SampleLevel(GSampler, float2(min(NoV, 0.999f), Roughness), 0.0f).rg;
		Specular = PrefilteredColor * (F * EnvBRDF.x + EnvBRDF.y);
	}
	else if (GPerFrameCB.IBLMode == IBL_MODE_SPHERICAL_HARMONICS_EXPONENTIAL)
	{
	}
	else if (GPerFrameCB.IBLMode == IBL_MODE_REFERENCE)
	{
		uint NumSamples = 128;
		float3 Irradiance = 0.0f;

		for (uint DiffuseSampleIdx = 0; DiffuseSampleIdx < NumSamples; ++DiffuseSampleIdx)
		{
			float2 Xi = HaltonSequence2D(DiffuseSampleIdx + GPerFrameCB.NumFrames * NumSamples);
			float3 L = CosineSampleHemisphere(Xi, N);

			float NoL = saturate(dot(N, L));
			if (NoL > 0.0f)
			{
				// Standard way of cosine importance sampling, no need to multiply by NoL.
				Irradiance += GEnvMap.SampleLevel(GSampler, L, 0).rgb;
			}
		}

		for (uint SpecularSampleIdx = 0; SpecularSampleIdx < NumSamples; ++SpecularSampleIdx)
		{
			float2 Xi = HaltonSequence2D(SpecularSampleIdx + GPerFrameCB.NumFrames * NumSamples);
			float3 H = ImportanceSampleGGX(Xi, Roughness, N);
			float3 L = normalize(2.0f * dot(V, H) * H - V);

			float NoL = saturate(dot(N, L));
			float NoH = saturate(dot(N, H));
			float VoH = saturate(dot(V, H));

			if (NoL > 0.0f)
			{
				float3 Li = GEnvMap.SampleLevel(GSampler, L, 0).rgb;
				float G = GeometrySmith(NoL, NoV, Roughness);
				float G_Vis = G * VoH / (NoH * NoV);
				Specular += Li * F * G_Vis * NoL;
			}
		}

		Diffuse = (Irradiance * Albedo) * (1.0f / NumSamples);
		Specular *= (1.0f / NumSamples);
	}

	float3 Ambient = 0.0f;
	if (GPerFrameCB.MaterialMode == MATERIAL_MODE_DIFFUSE_AND_SPECULAR)
	{
		Ambient = KD * Diffuse + Specular;
	}
	else if (GPerFrameCB.MaterialMode == MATERIAL_MODE_DIFFUSE_ONLY)
	{
		Ambient = KD * Diffuse;
	}
	else if (GPerFrameCB.MaterialMode == MATERIAL_MODE_SPECULAR_ONLY)
	{
		Ambient = Specular;
	}
	Ambient *= AO;

	float3 Color = Ambient + Lo;

	if (GPerFrameCB.IBLMode == IBL_MODE_REFERENCE && GPerFrameCB.NumFrames > 0)
	{
		// Progressive rendering for reference mode
		uint2 Dimensions;
		GAccumulationBuffer.GetDimensions(Dimensions.x, Dimensions.y);
		float3 AccumulatedColor = GAccumulationBuffer.SampleLevel(GSampler, InPosition.xy / Dimensions, 0).rgb;
		float3 AccumulatedLinearColor = pow(AccumulatedColor, 2.2f);
		float3 AccumulatedHDRColor = AccumulatedLinearColor / (1.0f - AccumulatedLinearColor);
		Color = (AccumulatedHDRColor * (GPerFrameCB.NumFrames - 1) + Color) / GPerFrameCB.NumFrames;
	}

	Color = Color / (Color + 1.0f);
	Color = pow(Color, 1.0f / 2.2f);
	OutColor = float4(Color, 1.0f);
}
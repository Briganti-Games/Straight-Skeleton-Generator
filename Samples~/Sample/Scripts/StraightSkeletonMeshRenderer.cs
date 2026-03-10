
using System.Collections.Generic;
using UnityEngine;
using Briganti.StraightSkeletonGeneration;
using Unity.Mathematics;
using UnityEngine.Serialization;

public class StraightSkeletonMeshRenderer : MonoBehaviour
{
	public float2[] polygonClockwise; // inspector
	public float2[] holeCounterClockwise; // inspector

	public MeshFilter meshFilter; // inspector


	public void Start()
	{
		// set up the input structure
		PolygonWithHoles polygonWithHoles = new PolygonWithHoles(polygonClockwise, new List<float2[]>() { holeCounterClockwise });

		// generate the straight skeleton
		StraightSkeletonGenerator generator = new StraightSkeletonGenerator(polygonWithHoles, float.MaxValue, false);
		StraightSkeleton skeleton = generator.Generate();

		// convert into a mesh and render it
		StraightSkeletonMeshGenerator meshGenerator = new StraightSkeletonMeshGenerator(skeleton);
		Mesh mesh = meshGenerator.GenerateRoofMesh();
		meshFilter.sharedMesh = mesh;
	}
}
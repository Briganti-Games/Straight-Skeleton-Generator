
using System.Collections.Generic;
using UnityEngine;
using Briganti.StraightSkeletonGeneration;
using Unity.Mathematics;

public class StraightSkeletonMeshRenderer : MonoBehaviour
{
	public float2[] polygonCounterClockwise; // inspector
	public float2[] holeClockwise; // inspector

	public MeshFilter meshFilter; // inspector


	public void Start()
	{
		// set up the input structure
		PolygonWithHoles polygonWithHoles = new PolygonWithHoles(polygonCounterClockwise, new List<float2[]>() { holeClockwise });

		// generate the straight skeleton
		StraightSkeletonGenerator generator = new StraightSkeletonGenerator(polygonWithHoles, float.MaxValue, false);
		StraightSkeleton skeleton = generator.Generate();

		// convert into a mesh and render it
		StraightSkeletonMeshGenerator meshGenerator = new StraightSkeletonMeshGenerator(skeleton);
		Mesh mesh = meshGenerator.GenerateRoofMesh();
		meshFilter.sharedMesh = mesh;
	}
}
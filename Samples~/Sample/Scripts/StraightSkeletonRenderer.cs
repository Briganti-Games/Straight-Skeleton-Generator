
using System.Collections.Generic;
using UnityEngine;
using Briganti.StraightSkeletonGeneration;
using Unity.Mathematics;

public class StraightSkeletonRenderer : MonoBehaviour
{
	public float2[] polygonCounterClockwise; // inspector
	public float2[] holeClockwise; // inspector

	public Transform polygonLinePrefab; // inspector
	public Transform straightSkeletonLinePrefab; // inspector


	public void Start()
	{
		// set up the input structure
		PolygonWithHoles polygonWithHoles = new PolygonWithHoles(polygonCounterClockwise, new List<float2[]>() { holeClockwise });

		// render the original polygon in 3D
		GeneratePolygonLines(polygonWithHoles);

		// generate the straight skeleton
		StraightSkeletonGenerator generator = new StraightSkeletonGenerator(polygonWithHoles, float.MaxValue, false);
		StraightSkeleton skeleton = generator.Generate();

		// render it as a 3D object, using the straight skeleton's collapse time as the 3D height so it resembles a roof
		GenerateStraightSkeletonLines(skeleton);
	}

	private void GeneratePolygonLines(PolygonWithHoles polygon)
	{
		GeneratePolygonLines(polygon.outerContourCounterClockwise);

		foreach (float2[] hole in polygon.innerContoursClockwise) {
			GeneratePolygonLines(hole);
		}
	}

	private void GeneratePolygonLines(float2[] polygon)
	{
		for (int i = 0; i < polygon.Length; ++i)
		{
			Vector2 p1 = polygon[i];
			Vector2 p2 = polygon[(i + 1) % polygon.Length];

			AddLine(new Vector3(p1.x, 0, p1.y), new Vector3(p2.x, 0, p2.y), polygonLinePrefab);
		}
	}

	private void GenerateStraightSkeletonLines(StraightSkeleton skeleton)
	{
		for (int i = 0; i < skeleton.nEdges; ++i)
		{
			Edge edge = skeleton.edges[i];
			int vertexIndex1 = edge.prevVertexIndex;
			int vertexIndex2 = edge.nextVertexIndex;

			float vertexTime1 = skeleton.vertexTimes[vertexIndex1];
			float vertexTime2 = skeleton.vertexTimes[vertexIndex2];
			float2 vertexPos1 = skeleton.vertices[vertexIndex1];
			float2 vertexPos2 = skeleton.vertices[vertexIndex2];

			Vector3 p1 = new Vector3(vertexPos1.x, vertexTime1, vertexPos1.y);
			Vector3 p2 = new Vector3(vertexPos2.x, vertexTime2, vertexPos2.y);

			AddLine(p1, p2, straightSkeletonLinePrefab);
		}
	}

	private void AddLine(Vector3 p1, Vector3 p2, Transform prefab)
	{
		Transform line = Instantiate(prefab);

		// position in the center
		line.transform.position = (p1 + p2) * 0.5f;

		// scale the proper length
		line.transform.localScale = new Vector3(line.transform.localScale.x, line.transform.localScale.y, (p2 - p1).magnitude);

		// calculate the quaternion to rotate from a line that goes from [0,0,0] to [0,0,1] to the actual line orientation
		Vector3 fromDirection = new Vector3(0, 0, 1);
		Vector3 toDorection = (p2 - p1).normalized;
		Quaternion rotation = Quaternion.FromToRotation(fromDirection, toDorection);
		line.transform.rotation = rotation;
	}
}
using System;
using Unity.Mathematics;
using UnityEngine;
using System.Linq;
#if UNITY_EDITOR
using UnityEditor;
#endif

namespace Briganti.StraightSkeletonGeneration
{
	public class StraightSkeleton : VertexGraph
	{
		public float[] vertexTimes { get; private set; }

		public float maxTime { get; private set; } = float.MinValue;

		public StraightSkeleton(int maxVertices, int maxEdges) : base(maxVertices, maxEdges)
		{
			vertexTimes = new float[maxVertices];
		}

		public int AddVertex(float2 pos, float time)
		{
			int vertexIndex = base.AddVertex(pos);
			vertexTimes[vertexIndex] = time;
			maxTime = Mathf.Max(maxTime, time);
			return vertexIndex;
		}

		protected override void IncreaseVertexCapacity(int nVertices)
		{
			base.IncreaseVertexCapacity(nVertices);

			int oldLength = vertexTimes.Length;

			float[] newVertexTimes = new float[nVertices];
			Array.Copy(vertexTimes, newVertexTimes, oldLength);
			vertexTimes = newVertexTimes;
		}

		public void DrawGizmos(Camera camera)
		{

			// calculate the distance from the camera of all the points, so we can apply a scale for some offsets based on the distance from the camera
			var points = Enumerable.Range(0, nVertices).Select(i =>
			{
				float2 vertex = vertices[i];
				return new Vector3(vertex.x, vertexTimes[i], vertex.y);
			});

			// based on the distance from the camera, we can now calculate an appropriate offset for all the labels
			float viewportOffset = 0.005f;
			float offset = WavefrontGraph.GetCameraDistanceIndepententOffset(camera, points, viewportOffset);

			Gizmos.color = Color.red;
			for (int i = 0; i < nEdges; ++i)
			{
				ref var edge = ref edges[i];

				float prevCreationTime = vertexTimes[edge.prevVertexIndex];
				float nextCreationTime = vertexTimes[edge.nextVertexIndex];

				var prevVertex = vertices[edge.prevVertexIndex];
				var nextVertex = vertices[edge.nextVertexIndex];

				Gizmos.DrawLine(new Vector3(prevVertex.x, prevCreationTime, prevVertex.y), new Vector3(nextVertex.x, nextCreationTime, nextVertex.y));

#if UNITY_EDITOR
				float2 midPoint = (prevVertex + nextVertex) * 0.5f;
				float2 normalBack = math.normalize(Geometry2D.Rotate90DegreesCounterClockwise(nextVertex - prevVertex));
				midPoint += normalBack * offset; // 0.1f
				float midCreationTime = (prevCreationTime + nextCreationTime) * 0.5f;
				Handles.Label(Gizmos.matrix.MultiplyPoint(new Vector3(midPoint.x, midCreationTime, midPoint.y)), "" + i, WavefrontGraph.edgeGuiStyle);
#endif
			}

			Gizmos.color = Color.cyan;
			for (int i = 0; i < nVertices; ++i)
			{
				float2 vertex = vertices[i];
				float creationTime = vertexTimes[i];

				Gizmos.DrawSphere(new Vector3(vertex.x, creationTime, vertex.y), 0.002f);

#if UNITY_EDITOR
				// render the number with a random direction offset
				var idPoint = vertex;
				const int nDivisions = 16;
				float angle = Mathf.PI * 2 * (i % nDivisions) / (float)nDivisions;
				float2 dir = new float2(math.cos(angle), math.sin(angle));
				idPoint = idPoint + dir * offset;
				Handles.Label(Gizmos.matrix.MultiplyPoint(new Vector3(idPoint.x, creationTime, idPoint.y)), "" + i, WavefrontGraph.vertexGuiStyle);
#endif
			}
		}
	}
}
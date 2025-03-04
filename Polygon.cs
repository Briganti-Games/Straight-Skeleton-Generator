using System;
using System.Buffers;
using System.Collections.Generic;

using Unity.Mathematics;

using UnityEngine;
using UnityEngine.Pool;

namespace Briganti.StraightSkeletonGeneration
{
	public class Polygon
	{
		public List<Vector2> vertices;

		public Polygon(in ReadOnlySpan<float2> vertices)
		{
			List<Vector2> processedVertices = new(vertices.Length);

			// don't allow for duplicates
			for (int i = 0; i < vertices.Length; i++)
			{
				Vector2 vertex = vertices[i];
				if (processedVertices.Count > 0 && (processedVertices[processedVertices.Count - 1] - vertex).sqrMagnitude < 0.0001f) continue;
				processedVertices.Add(vertex);
			}

			if (!IsCounterClockwise(processedVertices))
			{
				processedVertices.Reverse();
			}
			this.vertices = processedVertices;

			if (processedVertices.Count < 3) throw new ArgumentException($"A polygon must always consist of at least 3 vertices, but only these vertices remain after {string.Join(" -> ", processedVertices)}");
		}

		public Polygon(IEnumerable<Vector2> vertices)
		{
			List<Vector2> processedVertices = new();

			// don't allow for duplicates
			foreach (Vector2 vertex in vertices)
			{
				if (processedVertices.Count > 0 && (processedVertices[processedVertices.Count - 1] - vertex).sqrMagnitude < 0.0001f) continue;
				processedVertices.Add(vertex);
			}

			if (!IsCounterClockwise(processedVertices))
			{
				processedVertices.Reverse();
			}
			this.vertices = processedVertices;

			if (processedVertices.Count < 3) throw new ArgumentException($"A polygon must always consist of at least 3 vertices, but only these vertices remain after {string.Join(" -> ", processedVertices)}");
		}

		private static bool IsCounterClockwise(List<Vector2> points)
		{
			// https://www.element84.com/blog/determining-the-winding-of-a-polygon-given-as-a-set-of-ordered-points
			// https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
			float sum = 0.0f;
			for (int i = 0; i < points.Count; ++i)
			{
				sum += (points[(i + 1) % points.Count].x - points[i].x) * (points[(i + 1) % points.Count].y + points[i].y);
			}

			return sum < 0;
		}

		public bool IsConvex()
		{
			int sign = 0;
			for (int i = 0; i < vertices.Count - 1; ++i)
			{
				Vector2 p1 = vertices[i];
				Vector2 p2 = vertices[i + 1];
				Vector2 p3 = vertices[(i + 2) % vertices.Count];

				Vector2 v1 = p2 - p1;
				Vector2 v2 = p3 - p2;

				int newSign = (int)Mathf.Sign(Vector3.Cross(v1, v2).z);
				if (sign != 0 && sign != newSign) return false;
				sign = newSign;
			}

			return true;
		}

		public bool IsHorizontallyMonotone()
		{
			return IsMonotone(0);
		}

		public bool IsVerticallyMonotone()
		{
			return IsMonotone(1);
		}

		private bool IsMonotone(int dim)
		{
			// https://cs.stackexchange.com/questions/1577/how-do-i-test-if-a-polygon-is-monotone-with-respect-to-a-line

			bool IsLess(int i, int j)
			{
				return (vertices[i][dim] < vertices[j][dim]) || (vertices[i][dim] == vertices[j][dim] && i < j);
			}

			// see how many points are "leftmost" compared to their neighbours
			int localMinimums = 0;
			int n = vertices.Count;
			for (int i = 0; i < n; ++i)
			{
				if (IsLess(i, (i - 1 + n) % n) && IsLess(i, (i + 1) % n))
				{
					//Debug.Log("Vertex #" + i + " " + vertices[i] + " is a local minimum because the surrounding vertices are " + vertices[(i - 1 + n) % n] + " and " + vertices[(i + 1) % n]);
					++localMinimums;
					if (localMinimums >= 2)
					{
						return false;
					}
				}
			}

			return true;
		}

		public Mesh ToMesh(int[] indices)
		{
			Mesh mesh = new();

			int vertexCount = this.vertices.Count;
			int indexCount = indices.Length;
			var vertices = ArrayPool<Vector3>.Shared.Rent(vertexCount);
			var normals = ArrayPool<Vector3>.Shared.Rent(vertexCount);
			var triangles = ListPool<int>.Get();
			if (triangles.Capacity < indexCount)
				triangles.Capacity = indexCount;
			try
			{
				for (int i = 0; i < vertexCount; ++i)
				{
					Vector2 v = this.vertices[i];
					vertices[i] = new Vector3(v.x, 0, v.y);
					normals[i] = new Vector3(0, 1, 0);
				}

				// invert the indices to clockwise
				for (int i = 0; i < indexCount; i += 3)
				{
					triangles.Add(indices[i + 0]);
					triangles.Add(indices[i + 2]);
					triangles.Add(indices[i + 1]);
				}

				mesh.SetVertices(vertices, 0, vertexCount);
				mesh.SetNormals(normals, 0, vertexCount);
				mesh.SetIndices(triangles, MeshTopology.Triangles, 0);
			}
			finally
			{
				ArrayPool<Vector3>.Shared.Return(vertices);
				ArrayPool<Vector3>.Shared.Return(normals);
				ListPool<int>.Release(triangles);
			}
			return mesh;
		}

		public override string ToString()
		{
			return "[" + string.Join(" -> ", vertices) + "]";
		}
	}
}
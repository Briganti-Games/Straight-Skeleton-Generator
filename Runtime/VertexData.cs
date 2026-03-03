using Unity.Mathematics;

namespace Briganti.StraightSkeletonGeneration
{
	public enum SplitPoint
	{
		PrevVertex,
		Edge,
		NextVertex,
	}

	public struct VertexData
	{
		public float2 velocity;
		public WavefrontVertexType type;
		public float creationTime;

		public int prevVertexIndex;
		public int nextVertexIndex;

		public int prevEdgeIndex;
		public int nextEdgeIndex;

		public bool inWavefront;

		// this is set by the reflex vertex itself
		public SplitPoint splitPoint;
		public bool partOfSplitEvent;
		public float2 splitPos;
		public float splitTime;
		public int splitEdge;
		public int queueId;

		// this is set after assigning the vertex to an edge event
		//public int nextSplitReflexVertexIndex;
		

		public void UpdateConnections(int prevVertexIndex, int nextVertexIndex, int prevEdgeIndex, int nextEdgeIndex)
		{
			this.prevVertexIndex = prevVertexIndex;
			this.nextVertexIndex = nextVertexIndex;
			this.prevEdgeIndex = prevEdgeIndex;
			this.nextEdgeIndex = nextEdgeIndex;
			this.inWavefront = true;
		}

		public override string ToString()
		{
			if (!inWavefront) return "";
			if (partOfSplitEvent) return $"Split event at time {splitTime} and pos {splitPos} of edge {splitEdge}";
			return $"{type} vertex with velocity {velocity} created at {creationTime}, prev vertex {prevVertexIndex}, next vertex {nextVertexIndex}, still active {inWavefront}";
		}
	}
}

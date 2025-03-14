using Unity.Mathematics;

namespace Briganti.StraightSkeletonGeneration
{
	public enum EventType
	{
		None,
		Edge,
		Split,
		Nook,
		NotInWavefront,
	}

	public static class EventTypeExtensions
	{
		public static bool IsBatchEvent(this EventType eventType)
		{
			return eventType == EventType.Edge || eventType == EventType.Split;
		}

		public static bool IsValid(this EventType eventType)
		{
			return eventType != EventType.None && eventType != EventType.NotInWavefront;
		}
	}


	public struct Vertex
	{
		public float2 pos;

		public Vertex(float2 pos)
		{
			this.pos = pos;
		}

		public override string ToString()
		{
			return $"{pos}";
		}
	}
}
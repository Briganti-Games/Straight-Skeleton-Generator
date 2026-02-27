using Briganti.StraightSkeletons.Priority_Queue;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;

namespace Briganti.StraightSkeletonGeneration
{
	public class StringBuilderStraightSkeletonLogger : IStraightSkeletonLogger
	{
		private StringBuilder ss = new StringBuilder();


		public void Log(string message, int depth)
		{
			for (int i = 0; i < depth; ++i) ss.Append('\t');
			ss.AppendLine(message);
		}

		public override string ToString()
		{
			return ss.ToString();
		}
	}
}
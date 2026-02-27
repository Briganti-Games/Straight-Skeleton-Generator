This is a C# implementation of the Straight Skeleton algorithm for Unity.

Straight Skeleton is an algorithm that takes a polygon (with holes) as an input and spits out a bunch of lines that resemble a roof of a building.

It basically turns this:

![A polygon with a hole in it](https://github.com/Briganti-Games/Straight-Skeleton-Generator/blob/master/Readme%20Images/polygon-input.jpg?raw=true)

Into this:

![A straight skeleton generated from a polygon](https://github.com/Briganti-Games/Straight-Skeleton-Generator/blob/master/Readme%20Images/straight-skeleton-output.jpg?raw=true)

This is very useful for procedurally generating roof structures.

There are several other implementations for straight skeletons out there, but all the ones I tried had several problems:
- Sometimes they were unstable for more complex polygons.
- Sometimes they were WAY too slow for realtime procedural generation.
- Sometimes they used too many resources or memory allocations.
- Sometimes they don't support polygons with holes.

This implementation does all of that. It has been battle-tested by thousands of users over the course of a year, because it's the core of the roof generation in [Dungeon Alchemist](https://www.dungeonalchemist.com/), the procedural map-generation app for roleplaying games that I'm working on. It's robust, works with very complex polygons with hundreds of vertices, is very fast and does not waste resources.

# Getting Started

Install this package using the Unity Package Manager via git URL.

First, go to **Window** > **Package Manager**.

Then, navigate to the following submenu:

![Screenshot of Unity Package Manager](https://github.com/Briganti-Games/Straight-Skeleton-Generator/blob/master/Readme%20Images/package-manager.png?raw=true)

Then, copy & paste the URL to this git package: https://github.com/Briganti-Games/Straight-Skeleton-Generator.git and press enter. This will start the import process.

Once imported, you can use the package.

This has been tested on Unity 2022.3, but should probably also work on earlier versions, as long as the Unity Mathematics package is supported, as this is the only dependency.

# How to use it?

If you want to see the example that was used to generate the screenshots above, you can import it from the Package Manager by going to the Samples tab after importing the package, and importing the sample:

![Screenshot of the sample import process in the Package Manager](https://github.com/Briganti-Games/Straight-Skeleton-Generator/blob/master/Readme%20Images/sample-import.png?raw=true)

Using it in code is very simple. You give it your polygon data, and it outputs a straight skeleton for the polygon. Here's an example:

```
// counter-clockwise
float2[] polygonContour = new float2[] { new float2(0, 0), new float2(0, 5), new float2(5, 5), new float2(5, 0) };

// holes must be defined clockwise
float2[] hole = new float2[] { new float2(2, 2), new float2(3, 2), new float2(3, 3), new float2(2, 3) };

// set up the input data
PolygonWithHoles input = new PolygonWithHoles(polygonContour, new List<float2[]>() { hole });

// generate the straight skeleton
StraightSkeletonGenerator generator = new StraightSkeletonGenerator(input);
StraightSkeleton straightSkeleton = generator.Generate();

// read the output edges...
for (int i = 0; i < straightSkeleton.nEdges; ++i)
{
    Edge edge = straightSkeleton.edges[i];
    // do something with this data
}
```

# Implementation Details

Straight skeletons were first introduced in a 1995 paper by Aichholzer et al.:
https://www.jucs.org/jucs_1_12/a_novel_type_of/Aichholzer_O.pdf

In this first paper, they already realized the potential of this algorithm to be used as a procedural roof generation system.

I started implementing this algorithm using the motorcycle graph algorithm, which is the academic state of the art for this subject:
https://www.sthu.org/research/straightskeleton/

However, I realized this algorithm is complete overkill for what I needed and is incredibly complex to implement. Instead, I resorted to a more classical approach using a priority queue, loosely based on several other papers:
https://www.cgal.org/Events/UserWorkshop/2004/straight_skeleton.pdf
https://www.phlatforum.com/xenforo/converted_files/12703=3958-Straight%20Skeletons%20Implementation.pdf

None of these papers actually get into the nitty gritty of making this work for arbitrary polygons, and there are an incredible amount of very specific edge cases to iron out, as well as many issues with floating point precision. It took several months of development, as well as a year and thousands of users's worth of feedback to find most of these edge cases and bugs. But I feel like the current state of the implementation is very robust and stable, and works well even for very complex polygons.

One thing to note is that there are some arbitrary EPS values defined in Geometry2D.cs, which strongly depend on the scale of the polygons you are using. These values were fine-tuned to work well for polygons that work in a range of [0,100], but if you are using different scales (for example, vertices ranging from 0 to 10000) I would recommend adjusting these epsilon values accordingly. They are used in the code to determine whether lines intersect, and they can make a real difference in the result for very complex polygons.

//=============================================================================================================
/**
* @file     projectingkdtree.cpp
* @author   Lars Debor <lars.debor@tu-ilmenau.de>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     May, 2017
*
* @section  LICENSE
*
* Copyright (C) 2017, Lars Debor and Matti Hamalainen. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that
* the following conditions are met:
*     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
*       following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
*       the following disclaimer in the documentation and/or other materials provided with the distribution.
*     * Neither the name of MNE-CPP authors nor the names of its contributors may be used
*       to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
* PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*
* @brief    ProjectingKdTree class definition.
*
*/


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "projectingkdtree.h"
#include<mne/mne_bem_surface.h>

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <algorithm>
#include <cmath>
#include <iostream>

//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace GEOMETRYINFO;



//*************************************************************************************************************
//=============================================================================================================
// DEFINE GLOBAL METHODS
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================




ProjectingKdTree::ProjectingKdTree(const MNELIB::MNEBemSurface &inSurface, quint32 bucketSize) : m_surface(inSurface), m_maxBucketSize(bucketSize)
{
    m_vertIndices = std::vector<qint32>(inSurface.rr.rows());
    std::iota(m_vertIndices.begin(), m_vertIndices.end(), 0);

    mergeSort(m_vertIndices, 0, 200000, 0);

//    for (int i = 0; i <= 2000; i++)
//    {
//        std::cout << m_vertIndices[i] << ";" << m_surface.rr(m_vertIndices[i],0) << " ";
//    }
    std::cout << '\n';

    if(bucketSize > 0 )
    {
        m_root = recursiveBuild(m_vertIndices.begin(), m_vertIndices.end(), 0);
    }
    else
    {
        std::cout << "ERROR: BucketSize of projecting tree is = 0\n";
    }
}

ProjectingKdTree::~ProjectingKdTree()
{
    recursiveClear(m_root);
}

qint32 ProjectingKdTree::findNearestNeighbor(const Eigen::Vector3d &sensorPosition) const
{
    qint32 champion = -1;
    double minDistance = std::numeric_limits<double>::max();
    recursiveSearch(sensorPosition, m_root, champion, minDistance);
    if(champion < 0)
    {
        std::cout << "ERROR: No neighbor found!\n";
    }
    return champion;
}



//*************************************************************************************************************

ProjectingKdTree::ProjectingNode* ProjectingKdTree::recursiveBuild(std::vector<qint32>::iterator bucketBegin, std::vector<qint32>::iterator bucketEnd, qint32 depth)
{

    ProjectingNode *nodePtr = new ProjectingNode;
    //leaf reached
    const qint32 numPoints = std::distance(bucketBegin, bucketEnd);
    if(numPoints <= m_maxBucketSize)
    {

        nodePtr->m_bucketBegin = bucketBegin;
        nodePtr->m_bucketEnd = bucketEnd;
        //mark as leafnode
        nodePtr->m_subTrees[0] = nullptr;
        nodePtr->m_subTrees[1] = nullptr;

    }
    //nonleaf
    else
    {
        //create subtrees
        const qint8 axis = depth % 3;
        nodePtr->m_divAxis = axis;

        std::sort(bucketBegin, bucketEnd, [&](qint32 lhs, qint32 rhs)
        {
            return m_surface.rr(lhs, axis) < m_surface.rr(rhs, axis);
        });

        const qint32 mid = numPoints / 2;
        //Median of 3
        const double median = (m_surface.rr(*bucketBegin, axis) + m_surface.rr(*(bucketEnd - 1), axis) + m_surface.rr(*(bucketBegin + mid), axis)) / 3;

        auto midIt = std::upper_bound(bucketBegin, bucketEnd, median, [&](const double &value, const qint32 &indx){
            return value < m_surface.rr(indx, axis);
        });

        nodePtr->m_divValue = median;
        //left subtree
        nodePtr->m_subTrees[0] = recursiveBuild(bucketBegin, midIt, depth + 1);
        //right subtree
        nodePtr->m_subTrees[1] = recursiveBuild(midIt, bucketEnd, depth + 1);
    }
    return nodePtr;

    //this node is a leaf: bucket + no subtrees
//    if(depth == maxDepth)
//    {

//        nodePtr->m_bucketPtr = vertIndices;
//        nodePtr->m_bucketSize = numPoints;
//        nodePtr->m_subTrees[0] =
//        nodePtr->m_subTrees[1] = nullptr;
//    }
    //subtrees + no bucket
//    else
//    {

//        //TO DO  pivot element with random
//        const qint32 pivot = (numPoints - 1) / 2;

//        //partition
//        std::nth_element(vertIndices, vertIndices + pivot, vertIndices + numPoints, [&](qint32 lhs, qint32 rhs)
//        {
//            return m_surface.rr(lhs, axis) < m_surface.rr(rhs, axis);
//        });


//        nodePtr->m_vertIndex = vertIndices[pivot];

//        //left subtree
//        nodePtr->m_subTrees[0] = recursiveBuild(vertIndices, pivot, depth + 1, maxDepth);
//        //right subtree
//        nodePtr->m_subTrees[1] = recursiveBuild(vertIndices + pivot, numPoints - pivot - 1, depth + 1, maxDepth);
//    }

//    return nodePtr;

}

void ProjectingKdTree::recursiveSearch(const Eigen::Vector3d &sensorPosition, ProjectingNode *node, qint32 &champion, double &minDistance) const
{
    //const qint32 index = node->m_vertIndex;
    //Leaf reached?
    if(node->m_subTrees[0] == nullptr && node->m_subTrees[1] == nullptr)
    {
        //lin search inside the bucket
        std::for_each(node->m_bucketBegin, node->m_bucketEnd, [&](qint32 &n){
            const double dist = distance3D(sensorPosition, n);
            if(dist < minDistance)
            {
                champion = n;
                minDistance = dist;
            }
        });
        return;
//        //lin search in bucket
//        for(qint32 i = 0; i < node->m_bucketSize; ++i)
//        {
//            const double dist = distance3D(sensorPosition, *node->m_bucketPtr + i);
//            if(dist < minDistance)
//            {
//                champion = *node->m_bucketPtr + i;
//                minDistance = dist;
//            }
//        }

    }

    //NonLeaf
//    const double dist = distance3D(sensorPosition, );
//    if(dist < minDistance)
//    {
//        minDistance = dist;
//        champion = index;
//    }
    const qint8 axis = node->m_divAxis;
    const double divValue = node->m_divValue;
    // left or right first ?
    const bool subTreeIndx = sensorPosition[axis] <= divValue ? 0 : 1;
    recursiveSearch(sensorPosition, node->m_subTrees[subTreeIndx], champion, minDistance);

    const double diff = std::fabs(sensorPosition[axis] - divValue);
    //search the other subtree if needed
    if(diff < minDistance)
    {
        recursiveSearch(sensorPosition, node->m_subTrees[!subTreeIndx], champion, minDistance);
    }
//    else
//    {
//        static int count = 0;
//        std::cout << "b" << count++ << std::endl;
    //    }
}

void ProjectingKdTree::recursiveClear(ProjectingKdTree::ProjectingNode *node)
{
    if(node == nullptr)
    {
        return;
    }
    recursiveClear(node->m_subTrees[0]);
    recursiveClear(node->m_subTrees[1]);
    delete node;
}

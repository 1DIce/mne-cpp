//=============================================================================================================
/**
* @file     projectingkdtree.h
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
* @brief     ProjectingKdTree class declaration.
*
*/

#ifndef GEOMETRYINFO_PROJECTINGKDTREE_H
#define GEOMETRYINFO_PROJECTINGKDTREE_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "geometryinfo_global.h"
#include <mne/mne_bem_surface.h>
#include <vector>
#include <cmath>
#include <algorithm>

//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>
//#include <QVector>


//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================

#include <Eigen/core>

//*************************************************************************************************************
//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

//namespace MNELIB {
//    class MNEBemSurface;
//}

//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE GEOMETRYINFO
//=============================================================================================================

namespace GEOMETRYINFO {


//*************************************************************************************************************
//=============================================================================================================
// GEOMETRYINFO FORWARD DECLARATIONS
//=============================================================================================================


//=============================================================================================================
/**
* Description of what this class is intended to do (in detail).
*
* @brief Brief description of this class.
*/

class GEOMETRYINFOSHARED_EXPORT ProjectingKdTree
{

public:
    typedef QSharedPointer<ProjectingKdTree> SPtr;            /**< Shared pointer type for ProjectingKdTree. */
    typedef QSharedPointer<const ProjectingKdTree> ConstSPtr; /**< Const shared pointer type for ProjectingKdTree. */

    //=========================================================================================================
    /**
    * Constructs a ProjectingKdTree object.
    */
    ProjectingKdTree(const MNELIB::MNEBemSurface &inSurface, quint32 bucketSize);

    //custom destructor
    ~ProjectingKdTree();

    //copy operations
    ProjectingKdTree(const ProjectingKdTree&) = delete;
    ProjectingKdTree& operator=(ProjectingKdTree&) = delete;

    /**
     * @brief findNearestNeighbor
     * @param sensorPosition
     * @return index of the nearest neighbor
     */
    //return -1 if m_root == nullptr
    qint32 findNearestNeighbor(const Eigen::Vector3d &sensorPosition) const;


protected:

private:

    //declares the nodes the tree is made of
    struct ProjectingNode
    {

        qint32 m_bucketBegin;
        qint32 m_bucketEnd;

        qint8 m_divAxis;
        double m_divValue;

        ProjectingNode *m_subTrees[2];
    };
    inline double distance3D(const Eigen::Vector3d &sensorPosition, qint32 vertIndex) const;
    ProjectingNode *recursiveBuild(qint32 bucketBegin, qint32 bucketEnd, qint32 depth);
    void recursiveSearch(const Eigen::Vector3d &sensorPosition, ProjectingNode *node, qint32 &champion, double &minDistance) const;
    void recursiveClear(ProjectingNode *node);
    //saves a pointer to the root node of the search tree
    ProjectingNode *m_root;
    const MNELIB::MNEBemSurface &m_surface;
    //all indices of the MENBemSurface rr are stored here
    std::vector<qint32> m_vertIndices;
    std::vector<std::vector<qint32>> m_sortedIndices;
    std::vector<qint32> m_medianVec;
    qint32 m_maxBucketSize;

    inline void mergeSort(std::vector<qint32> &m_vertIndices, qint32 start, qint32 end, qint8 axis);
    inline void merge (std::vector<qint32> &m_vertIndices, qint32 start, qint32 middle, qint32 end, qint8 axis);
    inline qint32 medianSort (std::vector<std::vector<qint32>> &vertIndices, qint32 start, qint32 end, qint8 axis);
    ProjectingNode *balancedTreeBuild(std::vector<qint32>::iterator bucketBegin, std::vector<qint32>::iterator bucketEnd, qint32 depth);
    inline  void stableSortAxis(std::vector<qint32>::iterator begin, std::vector<qint32>::iterator end, qint8 axis);
};


//*************************************************************************************************************
//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================
inline double ProjectingKdTree::distance3D(const Eigen::Vector3d &sensorPosition, qint32 vertIndex) const
{
    return sqrt(pow(m_surface.rr(vertIndex, 0) - sensorPosition[0], 2)  // x-cord
            + pow(m_surface.rr(vertIndex, 1) - sensorPosition[1], 2)    // y-cord
            + pow(m_surface.rr(vertIndex, 2) - sensorPosition[2], 2));  // z-cord
}

inline qint32 ProjectingKdTree::medianSort (std::vector<std::vector<qint32>> &vertIndices, qint32 start, qint32 end, qint8 axis)
{
    //copy vertIndices X -> Temp

    //nur teil array
    std::vector<qint32> tempVec = vertIndices[0];

    const qint32 mid = ceil((start+end)/2);
    const qint32 medianIdx = vertIndices[0][mid];
    //qint32 m_start = start;
    qint32 med = mid + 1;
    qint32 begin = start;
    qint32 otherAxis = (axis + 1) % 3;


    for (int i = start; i <= end; i++)
    {
       if (vertIndices[1][i] == medianIdx)
       {
           continue;
       }
       else
       {
           qint32 tempAxis = axis;
           while(m_surface.rr(vertIndices[1][i],tempAxis) == m_surface.rr(medianIdx,tempAxis))
           {
               tempAxis++;
               tempAxis = tempAxis % 3;
           }
           if (m_surface.rr(vertIndices[1][i],tempAxis) < m_surface.rr(medianIdx,tempAxis))
           {
               vertIndices[0][begin] = vertIndices[1][i];
               begin++;
           }
           else
           {
               vertIndices[0][med] = vertIndices[1][i];
               med++;
           }
       }

    }
    /////////////////////
    med = mid + 1;
    begin = start;
    //axis = (axis + 1) % 3;
    otherAxis = (axis + 1) % 3;

    for (int i = start; i <= end; i++)
    {
       if (vertIndices[2][i] == medianIdx)
       {
           continue;
       }
       else
       {
           qint32 tempAxis = axis;
           while(m_surface.rr(vertIndices[2][i],tempAxis) == m_surface.rr(medianIdx,tempAxis))
           {
               tempAxis++;
               tempAxis = tempAxis % 3;
           }
           if (m_surface.rr(vertIndices[2][i],tempAxis) < m_surface.rr(medianIdx,tempAxis))
           {
               vertIndices[1][begin] = vertIndices[2][i];
               begin++;
           }
           else
           {
               vertIndices[1][med] = vertIndices[2][i];
               med++;
           }
       }

    }
    //todo tempvec zurück
    vertIndices[2]= std::move(tempVec);
    return mid;

//    switch(axis)
//    {
//        case 0: qint32 median = m_vertIndicesX[ceil((start+end)/2)];
//            break;
//        case 1: qint32 median = m_vertIndicesY[ceil((start+end)/2)];
//            break;
//        case 2: qint32 median = m_vertIndicesZ[ceil((start+end)/2)];
//            break;
//        default: std::cout << "ERROR";
//    }

    //copy vertIndices Temp -> Z
}

inline void ProjectingKdTree::stableSortAxis(std::vector<qint32>::iterator begin, std::vector<qint32>::iterator end, qint8 axis)
{
    std::stable_sort(begin, end, [&](qint32 lhs, qint32 rhs){
        qint8 tempAxis = axis;
        while(m_surface.rr(lhs, tempAxis) == m_surface.rr(rhs, tempAxis))
        {
            tempAxis++;
            tempAxis = tempAxis % 3;
        }
        return m_surface.rr(lhs, tempAxis) < m_surface.rr(rhs, tempAxis);
    });
}

inline void ProjectingKdTree::mergeSort(std::vector<qint32> &m_vertIndices, qint32 start, qint32 end, qint8 axis)
{
    //std::cout << start << " + " << end << '\n' ;
    if (start < end)
    {
        int middle = int((end + start) / 2);
        mergeSort(m_vertIndices, start, middle, axis);
        mergeSort(m_vertIndices, middle+1, end, axis);
        merge(m_vertIndices, start, middle, end, axis);
    }
    return;
}

inline void ProjectingKdTree::merge (std::vector<qint32> &m_vertIndices, qint32 start, qint32 middle, qint32 end, qint8 axis)
{
//    std::vector<qint32> m_left (middle-start+1);
//    for (int l_counter = 0; l_counter < middle-start; l_counter++)
//    {
//        m_left[l_counter] = m_vertIndices[start+l_counter];
//        std::cout << m_left[l_counter] << " ";
//    }

//    m_left.push_back(100000);

//    std::cout << '\n';
//    std::vector<qint32> m_right (end-middle+2);
//    for (int r_counter = 0; r_counter < end-middle+1; r_counter++)
//    {
//        m_right[r_counter] = m_vertIndices[middle+r_counter];
//        std::cout << m_right[r_counter] << " ";
//    }

//    m_right.push_back(100000);

//    int m_cmp_i = 0;
//    int m_cmp_j = 0;

//    for (int k = start; k <= end; k++)
//    {
//        if (m_surface.rr(m_left[m_cmp_i],0) == m_surface.rr(m_right[m_cmp_j],0))
//        {
//            if(m_surface.rr(m_left[m_cmp_i],1) == m_surface.rr(m_right[m_cmp_j],1))
//            {
//                if (m_surface.rr(m_left[m_cmp_i],2) < m_surface.rr(m_right[m_cmp_j],2))
//                {
//                    m_vertIndices[k] = m_left[m_cmp_i];
//                    m_cmp_i++;
//                }
//                else
//                {
//                    m_vertIndices[k] = m_right[m_cmp_j];
//                    m_cmp_j++;
//                }
//            }
//            else
//            {
//                if (m_surface.rr(m_left[m_cmp_i],1) < m_surface.rr(m_right[m_cmp_j],1))
//                {
//                    m_vertIndices[k] = m_left[m_cmp_i];
//                    m_cmp_i++;
//                }
//                else
//                {
//                    m_vertIndices[k] = m_right[m_cmp_j];
//                    m_cmp_j++;
//                }
//            }
//        }
//        else
//        {
//            if (m_surface.rr(m_left[m_cmp_i],0) < m_surface.rr(m_right[m_cmp_j],0))
//            {
//                m_vertIndices[k] = m_left[m_cmp_i];
//                m_cmp_i++;
//            }
//            else
//            {
//                m_vertIndices[k] = m_right[m_cmp_j];
//                m_cmp_j++;
//            }

//        }
//    }

    std::vector<qint32> m_temp (200000);
    qint32 m_i = start;
    qint32 m_k = start;
    qint32 m_j = middle+1;

    switch (axis)
    {
        case 0:
        while (m_i<= middle && m_j <= end)
                        {
                        if (m_surface.rr(m_vertIndices[m_i],0) == m_surface.rr(m_vertIndices[m_j],0))
                        {
                            if (m_surface.rr(m_vertIndices[m_i],1) == m_surface.rr(m_vertIndices[m_j],1))
                            {
                                if (m_surface.rr(m_vertIndices[m_i],2) < m_surface.rr(m_vertIndices[m_j],2))
                                {
                                    m_temp[m_k] = m_vertIndices[m_i];
                                    m_k++;
                                    m_i++;
                                }
                                else
                                {
                                    m_temp[m_k] = m_vertIndices[m_j];
                                    m_k++;
                                    m_j++;
                                }
                            }
                            else
                            {
                                if (m_surface.rr(m_vertIndices[m_i],1) < m_surface.rr(m_vertIndices[m_j],1))
                                {
                                    m_temp[m_k] = m_vertIndices[m_i];
                                    m_k++;
                                    m_i++;
                                }
                                else
                                {
                                    m_temp[m_k] = m_vertIndices[m_j];
                                    m_k++;
                                    m_j++;
                                }
                            }
                        }
                        else
                        {
                            if (m_surface.rr(m_vertIndices[m_i],0) < m_surface.rr(m_vertIndices[m_j],0))
                            {
                                m_temp[m_k] = m_vertIndices[m_i];
                                m_k++;
                                m_i++;
                            }
                            else
                            {
                                m_temp[m_k] = m_vertIndices[m_j];
                                m_k++;
                                m_j++;
                            }
                        }
                    }

                    while (m_i <= middle)
                    {
                        m_temp[m_k] = m_vertIndices[m_i];
                        m_k++;
                        m_i++;
                    }

                    while (m_j <= end)
                    {
                        m_temp[m_k] = m_vertIndices[m_j];
                        m_k++;
                        m_j++;
                    }

                    for (int i = start; i < m_k; i++)
                    {
                        m_vertIndices[i] = m_temp[i];
                    }
            break;
        case 1:
        while (m_i<= middle && m_j <= end)
                        {
                        if (m_surface.rr(m_vertIndices[m_i],1) == m_surface.rr(m_vertIndices[m_j],1))
                        {
                            if (m_surface.rr(m_vertIndices[m_i],2) == m_surface.rr(m_vertIndices[m_j],2))
                            {
                                if (m_surface.rr(m_vertIndices[m_i],0) < m_surface.rr(m_vertIndices[m_j],0))
                                {
                                    m_temp[m_k] = m_vertIndices[m_i];
                                    m_k++;
                                    m_i++;
                                }
                                else
                                {
                                    m_temp[m_k] = m_vertIndices[m_j];
                                    m_k++;
                                    m_j++;
                                }
                            }
                            else
                            {
                                if (m_surface.rr(m_vertIndices[m_i],2) < m_surface.rr(m_vertIndices[m_j],2))
                                {
                                    m_temp[m_k] = m_vertIndices[m_i];
                                    m_k++;
                                    m_i++;
                                }
                                else
                                {
                                    m_temp[m_k] = m_vertIndices[m_j];
                                    m_k++;
                                    m_j++;
                                }
                            }
                        }
                        else
                        {
                            if (m_surface.rr(m_vertIndices[m_i],1) < m_surface.rr(m_vertIndices[m_j],1))
                            {
                                m_temp[m_k] = m_vertIndices[m_i];
                                m_k++;
                                m_i++;
                            }
                            else
                            {
                                m_temp[m_k] = m_vertIndices[m_j];
                                m_k++;
                                m_j++;
                            }
                        }
                    }

                    while (m_i <= middle)
                    {
                        m_temp[m_k] = m_vertIndices[m_i];
                        m_k++;
                        m_i++;
                    }

                    while (m_j <= end)
                    {
                        m_temp[m_k] = m_vertIndices[m_j];
                        m_k++;
                        m_j++;
                    }

                    for (int i = start; i < m_k; i++)
                    {
                        m_vertIndices[i] = m_temp[i];
                    }
            break;
        case 2:
        while (m_i<= middle && m_j <= end)
                        {
                        if (m_surface.rr(m_vertIndices[m_i],2) == m_surface.rr(m_vertIndices[m_j],2))
                        {
                            if (m_surface.rr(m_vertIndices[m_i],0) == m_surface.rr(m_vertIndices[m_j],0))
                            {
                                if (m_surface.rr(m_vertIndices[m_i],1) < m_surface.rr(m_vertIndices[m_j],1))
                                {
                                    m_temp[m_k] = m_vertIndices[m_i];
                                    m_k++;
                                    m_i++;
                                }
                                else
                                {
                                    m_temp[m_k] = m_vertIndices[m_j];
                                    m_k++;
                                    m_j++;
                                }
                            }
                            else
                            {
                                if (m_surface.rr(m_vertIndices[m_i],0) < m_surface.rr(m_vertIndices[m_j],0))
                                {
                                    m_temp[m_k] = m_vertIndices[m_i];
                                    m_k++;
                                    m_i++;
                                }
                                else
                                {
                                    m_temp[m_k] = m_vertIndices[m_j];
                                    m_k++;
                                    m_j++;
                                }
                            }
                        }
                        else
                        {
                            if (m_surface.rr(m_vertIndices[m_i],2) < m_surface.rr(m_vertIndices[m_j],2))
                            {
                                m_temp[m_k] = m_vertIndices[m_i];
                                m_k++;
                                m_i++;
                            }
                            else
                            {
                                m_temp[m_k] = m_vertIndices[m_j];
                                m_k++;
                                m_j++;
                            }
                        }
                    }

                    while (m_i <= middle)
                    {
                        m_temp[m_k] = m_vertIndices[m_i];
                        m_k++;
                        m_i++;
                    }

                    while (m_j <= end)
                    {
                        m_temp[m_k] = m_vertIndices[m_j];
                        m_k++;
                        m_j++;
                    }

                    for (int i = start; i < m_k; i++)
                    {
                        m_vertIndices[i] = m_temp[i];
                    }
            break;
        default: std::cout << "ERROR" ;
    }
}

} // namespace GEOMETRYINFO

#endif // GEOMETRYINFO_PROJECTINGKDTREE_H

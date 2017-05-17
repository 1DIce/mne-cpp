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

        std::vector<qint32>::iterator m_bucketBegin;
        std::vector<qint32>::iterator m_bucketEnd;

        qint8 m_divAxis;
        double m_divValue;

        ProjectingNode *m_subTrees[2];
    };
    inline double distance3D(const Eigen::Vector3d &sensorPosition, qint32 vertIndex) const;
    ProjectingNode *recursiveBuild(std::vector<qint32>::iterator bucketBegin, std::vector<qint32>::iterator bucketEnd, qint32 depth);
    void recursiveSearch(const Eigen::Vector3d &sensorPosition, ProjectingNode *node, qint32 &champion, double &minDistance) const;
    void recursiveClear(ProjectingNode *node);
    //saves a pointer to the root node of the search tree
    ProjectingNode *m_root;
    const MNELIB::MNEBemSurface &m_surface;
    //all indices of the MENBemSurface rr are stored here
    std::vector<qint32> m_vertIndices;
    qint32 m_maxBucketSize;

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

} // namespace GEOMETRYINFO

#endif // GEOMETRYINFO_PROJECTINGKDTREE_H

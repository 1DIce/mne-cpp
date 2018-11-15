//=============================================================================================================
/**
* @file     fiffrawmodel.cpp
* @author   Simon Heinke <simon.heinke@tu-ilmenau.de>;
*           Lars Debor <lars.debor@tu-ilmenau.de>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     October, 2018
*
* @section  LICENSE
*
* Copyright (C) 2018, Simon Heinke, Lars Debor and Matti Hamalainen. All rights reserved.
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
* @brief    FiffRawModel class definition.
*
*/


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "fiffrawmodel.h"
#include "../Utils/metatypes.h"

#include <algorithm>


//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtConcurrent/QtConcurrent>


//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace ANSHAREDLIB;
using namespace FIFFLIB;


//*************************************************************************************************************
//=============================================================================================================
// DEFINE GLOBAL METHODS
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

FiffRawModel::FiffRawModel(QObject *pParent)
    : AbstractModel(pParent)
{

}


//*************************************************************************************************************

FiffRawModel::FiffRawModel(QFile& inFile,
                           qint32 iSamplesPerBlock,
                           qint32 iWindowSize,
                           qint32 iPreloadBufferSize,
                           QObject *pParent)
    : AbstractModel(pParent),
      m_iSamplesPerBlock(iSamplesPerBlock),
      m_iWindowSize(iWindowSize),
      m_iPreloadBufferSize(iPreloadBufferSize),
      m_iTotalBlockCount(m_iWindowSize + 2 * m_iPreloadBufferSize),
      m_iFiffCursorBegin(-1),
      m_blockLoadFutureWatcher(),
      m_bCurrentlyLoading(false),
      m_dataMutex(),
      m_pFiffIO(),
      m_pFiffInfo(),
      m_ChannelInfoList()
{
    // connect data reloading: this is done concurrently
    connect(&m_blockLoadFutureWatcher,
            &QFutureWatcher<int>::finished,
            [this]() {
                postBlockLoad(m_blockLoadFutureWatcher.future().result());
            });


    initFiffData(inFile);
}


//*************************************************************************************************************

FiffRawModel::~FiffRawModel()
{

}


//*************************************************************************************************************

void FiffRawModel::initFiffData(QFile &inFile)
{
    // build FiffIO
    m_pFiffIO = QSharedPointer<FiffIO>::create(inFile);

    // load channel infos
    for(qint32 i=0; i < m_pFiffIO->m_qlistRaw[0]->info.nchan; ++i)
        m_ChannelInfoList.append(m_pFiffIO->m_qlistRaw[0]->info.chs[i]);

    // load FiffInfo
    m_pFiffInfo = FiffInfo::SPtr(new FiffInfo(m_pFiffIO->m_qlistRaw[0]->info));

    // build datastructure that is to be filled with data from the file
    MatrixXd data, times;
    // append a matrix pair for each block
    for(int i = 0; i < m_iTotalBlockCount; ++i) {
        m_lData.append(QSharedPointer<QPair<MatrixXd, MatrixXd> >::create(qMakePair(data, times)));
    }

    if(m_pFiffIO->m_qlistRaw.empty()) {
        qDebug() << "[FiffRawModel::loadFiffData] File " << inFile.fileName() << " does not contain any Fiff data";
        return;
    }

    // Fiff file is not empty, set cursor somewhere into Fiff file
    m_iFiffCursorBegin = m_pFiffIO->m_qlistRaw[0]->first_samp;

    int start = m_iFiffCursorBegin;
    // for some reason the read_raw_segment function works with inclusive upper bound
    int end = start + m_iSamplesPerBlock - 1;

    // read in all blocks, use the already prepaired list m_lData
    for(auto &pairPointer : m_lData) {
        if(m_pFiffIO->m_qlistRaw[0]->read_raw_segment(pairPointer->first, pairPointer->second, start, end)) {
            // qDebug() << "[FiffRawmodel::loadFiffData] Successfully read a block ";
        } else {
            qDebug() << "[FiffRawModel::loadFiffData] Could not read samples " << start << " to " << end;
            return;
        }

        start += m_iSamplesPerBlock;
        end += m_iSamplesPerBlock;
    }

    qDebug() << "[FiffRawModel::initFiffData] Loaded " << m_lData.size() << " blocks";
}


//*************************************************************************************************************

QVariant FiffRawModel::data(const QModelIndex &index, int role) const
{
    // early filtering of unimplemented display roles
    if (role != Qt::DisplayRole) {
        qDebug() << "[FiffRawModel] Role " << role << " not implemented yet !";
        return QVariant();
    }

    if (index.isValid()) {
        // channel names
        if(index.column() == 0) {
            return QVariant(m_ChannelInfoList[index.row()].ch_name);
        }

        // data
        if (index.column() == 1) {
            QVariant result;

            switch (role) {
            case Qt::DisplayRole:
                // in order to avoid extensive copying of data, we take advantage of Eigen matrices being organized row-wise
                // see ChannelData.h for more details

                // wait until its save to access data (that is if no data insertion is going on right now)
                m_dataMutex.lock();

                // wrap in ChannelData container and then wrap into QVariant
                result.setValue(ChannelData(m_lData, index.row()));

                m_dataMutex.unlock();

                return result;
                break;
            }
        }
        if (index.column() != 0 && index.column() != 1) {
            qDebug() << "[FiffRawModel] Column " << index.column() << " not implemented !";
            return QVariant();
        }
    }

    qDebug() << "[FiffRawModel::data] Warning, non of the presumed cases took effect";
    return QVariant();
}


//*************************************************************************************************************

Qt::ItemFlags FiffRawModel::flags(const QModelIndex &index) const
{
    // TODO implement stuff
    return QAbstractItemModel::flags(index);
}


//*************************************************************************************************************

QModelIndex FiffRawModel::index(int row, int column, const QModelIndex &parent) const
{
    return createIndex(row, column);
}


//*************************************************************************************************************

QModelIndex FiffRawModel::parent(const QModelIndex &index) const
{
    // TODO implement stuff
    return QModelIndex();
}


//*************************************************************************************************************

int FiffRawModel::rowCount(const QModelIndex &parent) const
{
    if(m_ChannelInfoList.empty() == false)
        return m_ChannelInfoList.size();
    else
        return 0;
}


//*************************************************************************************************************

int FiffRawModel::columnCount(const QModelIndex &parent) const
{
    // TODO implement stuff
    return 0;
}


//*************************************************************************************************************

bool FiffRawModel::hasChildren(const QModelIndex &parent) const
{
    // TODO implement stuff
    return false;
}


//*************************************************************************************************************

void FiffRawModel::updateScrollPosition(qint32 relativeFiffCursor)
{
    // check if we are currently loading something in the background. This is a rudimentary solution.
    if (m_bCurrentlyLoading) {
        // qDebug() << "[FiffRawModel::updateScrollPosition] Background operation still pending, try again later...";
        return;
    }

    qint32 targetCursor = absoluteFirstSample() + relativeFiffCursor;

    if (targetCursor < m_iFiffCursorBegin + m_iPreloadBufferSize * m_iSamplesPerBlock) {
        // time to move the loaded window. Calculate distance in blocks
        qint32 sampleDist = (m_iFiffCursorBegin + m_iPreloadBufferSize * m_iSamplesPerBlock) - targetCursor;
        qint32 blockDist = (qint32) ceil(((double) sampleDist) / ((double) m_iSamplesPerBlock));

        if (blockDist >= m_iTotalBlockCount) {
            // we must "jump" to the new cursor ...
            m_iFiffCursorBegin = std::max(absoluteFirstSample(), m_iFiffCursorBegin -  blockDist * m_iSamplesPerBlock);
            // ... and load the whole model anew
            startBackgroundOperation(&FiffRawModel::loadLaterBlocks, m_iTotalBlockCount);
        } else {
            // there are some blocks in the intersection of the old and the new window that can stay in the buffer:
            // simply load earlier blocks
            startBackgroundOperation(&FiffRawModel::loadEarlierBlocks, blockDist);
        }
    }
    else if (targetCursor >= m_iFiffCursorBegin + (m_iPreloadBufferSize + m_iWindowSize) * m_iSamplesPerBlock) {
        // time to move the loaded window. Calculate distance in blocks
        qint32 sampleDist = targetCursor - (m_iFiffCursorBegin + (m_iPreloadBufferSize + m_iWindowSize) * m_iSamplesPerBlock);
        qint32 blockDist = (qint32) ceil(((double) sampleDist) / ((double) m_iSamplesPerBlock));

        if (blockDist >= m_iTotalBlockCount) {
            // we must "jump" to the new cursor ...
            m_iFiffCursorBegin = std::min(absoluteLastSample() - m_iTotalBlockCount * m_iSamplesPerBlock, m_iFiffCursorBegin +  blockDist * m_iSamplesPerBlock);
            // ... and load the whole model anew
            startBackgroundOperation(&FiffRawModel::loadLaterBlocks, m_iTotalBlockCount);
        } else {
            // there are some blocks in the intersection of the old and the new window that can stay in the buffer:
            // simply load later blocks
            startBackgroundOperation(&FiffRawModel::loadLaterBlocks, blockDist);
        }
    }
}


//*************************************************************************************************************

void FiffRawModel::startBackgroundOperation(int (FiffRawModel::*loadFunction)(int), int iBlocksToLoad)
{
    m_bCurrentlyLoading = true;
    QFuture<int> future = QtConcurrent::run(this, loadFunction, iBlocksToLoad);
    m_blockLoadFutureWatcher.setFuture(future);
}


//*************************************************************************************************************

int FiffRawModel::loadEarlierBlocks(qint32 numBlocks)
{
    // check if start of file was reached:
    int leftSamples = (m_iFiffCursorBegin - numBlocks * m_iSamplesPerBlock) - absoluteFirstSample();
    if (leftSamples < 0) {
        qDebug() << "Reached start of file !";
        // see how many blocks we still can load
        int maxNumBlocks = (m_iFiffCursorBegin - absoluteFirstSample()) / m_iSamplesPerBlock;
        qDebug() << "Loading " << maxNumBlocks << " earlier blocks instead of requested " << numBlocks;
        if (maxNumBlocks != 0) {
            numBlocks = maxNumBlocks;
        }
        else {
            // nothing to be done, cant load any more blocks
            // return 0, meaning that this was a loading of earlier blocks
            return 0;
        }
    }
    // we expect m_lNewData to be empty:
    if (m_lNewData.isEmpty() == false) {
        qDebug() << "[FiffRawModel::loadEarlierBlocks] FATAL, temporary data storage non empty !";
        return -1;
    }
    // build data structures to be filled from file
    MatrixXd data, times;
    for(int i = 0; i < numBlocks; ++i) {
        m_lNewData.append(QSharedPointer<QPair<MatrixXd, MatrixXd> >::create(qMakePair(data, times)));
    }
    // initialize start and end indices
    int start = m_iFiffCursorBegin;
    int end = start + m_iSamplesPerBlock - 1;

    // read data, use the already prepared list m_lNewData
    for(auto &pairPointer : m_lNewData) {
        start -= m_iSamplesPerBlock;
        end -= m_iSamplesPerBlock;
        if(m_pFiffIO->m_qlistRaw[0]->read_raw_segment(pairPointer->first, pairPointer->second, start, end)) {
            // qDebug() << "[FiffRawModel::loadFiffData] Successfully read a block ";
        } else {
            qDebug() << "[FiffRawModel::loadEarlierBlocks] Could not read block ";
            return -1;
        }
    }

    // adjust fiff cursor
    m_iFiffCursorBegin = start;

    // return 0, meaning that this was a loading of earlier blocks
    return 0;
}


//*************************************************************************************************************

int FiffRawModel::loadLaterBlocks(qint32 numBlocks)
{
    // check if end of file is reached:
    int leftSamples = absoluteLastSample() - (m_iFiffCursorBegin + (m_iTotalBlockCount + numBlocks) * m_iSamplesPerBlock);
    if (leftSamples < 0) {
        qDebug() << "Reached end of file !";
        // see how many blocks we still can load
        int maxNumBlocks = (absoluteLastSample() - (m_iFiffCursorBegin + m_iTotalBlockCount * m_iSamplesPerBlock)) / m_iSamplesPerBlock;
        qDebug() << "Loading " << maxNumBlocks << " later blocks instead of requested " << numBlocks;
        if (maxNumBlocks != 0) {
            numBlocks = maxNumBlocks;
        }
        else {
            // nothing to be done, cant load any more blocks
            // return 1, meaning that this was a loading of later blocks
            return 1;
        }

    }
    // we expect m_lNewData to be empty:
    if (m_lNewData.isEmpty() == false) {
        qDebug() << "[FiffRawModel::loadLaterBlocks] FATAL, temporary data storage non empty !";
        return -1;
    }
    // build data structures to be filled from file
    MatrixXd data, times;
    for(int i = 0; i < numBlocks; ++i) {
        m_lNewData.append(QSharedPointer<QPair<MatrixXd, MatrixXd> >::create(qMakePair(data, times)));
    }
    // initialize start and end indices
    int start = m_iFiffCursorBegin + m_iTotalBlockCount * m_iSamplesPerBlock;
    int end = start + m_iSamplesPerBlock - 1;

    // read data, use the already prepaired list m_lNewData
    for(auto &pairPointer : m_lNewData) {
        if(m_pFiffIO->m_qlistRaw[0]->read_raw_segment(pairPointer->first, pairPointer->second, start, end)) {
            // qDebug() << "[FiffRawModel::loadFiffData] Successfully read a block ";
        } else {
            qDebug() << "[FiffRawModel::loadLaterBlocks] Could not read block ";
            return -1;
        }

        start += m_iSamplesPerBlock;
        end += m_iSamplesPerBlock;
    }

    // adjust fiff cursor
    m_iFiffCursorBegin += numBlocks * m_iSamplesPerBlock;

    // return 1, meaning that this was a loading of later blocks
    return 1;
}


//*************************************************************************************************************

void FiffRawModel::postBlockLoad(int result)
{
    switch(result){
    case -1:
        qDebug() << "[FiffRawModel::postBlockLoad] QFuture returned an error: " << result;
        break;
    case 0:
    {
        // insertion of earlier blocks
        int iNewBlocks = m_lNewData.size();

        m_dataMutex.lock();
        for (int i = 0; i < iNewBlocks; ++i) {
            m_lData.prepend(m_lNewData.first());
            m_lNewData.removeFirst();
            // @TODO check if this really frees the associated memory
            m_lData.removeLast();
        }
        m_dataMutex.unlock();

        break;
    }
    case 1:
    {
        // insertion of later blocks
        int iNewBlocks = m_lNewData.size();

        m_dataMutex.lock();
        for (int i = 0; i < iNewBlocks; ++i) {
            m_lData.append(m_lNewData.first());
            m_lNewData.removeFirst();
            m_lData.removeFirst();
        }
        m_dataMutex.unlock();

        break;
    }
    default:
        qDebug() << "[FiffRawModel::postBlockLoad] FATAL Non-intended return value: " << result;
    }

    m_bCurrentlyLoading = false;
}

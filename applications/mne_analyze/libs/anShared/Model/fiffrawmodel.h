//=============================================================================================================
/**
* @file     fiffrawmodel.h
* @author   Simon Heinke <simon.heinke@tu-ilmenau.de>;
*           Lars Debor <lars.debor@tu-ilmenau.de>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     April, 2018
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
* @brief     FiffRawModel class declaration.
*
*/

#ifndef ANSHAREDLIB_FIFFRAWMODEL_H
#define ANSHAREDLIB_FIFFRAWMODEL_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../anshared_global.h"
#include "../Utils/types.h"
#include "abstractmodel.h"

#include <fiff/fiff.h>
#include <fiff/fiff_io.h>

//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>
#include <QDebug>
#include <QFile>
#include <QFuture>
#include <QFutureWatcher>


//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE ANSHAREDLIB
//=============================================================================================================

namespace ANSHAREDLIB {


//*************************************************************************************************************
//=============================================================================================================
// ANSHAREDLIB FORWARD DECLARATIONS
//=============================================================================================================


//=============================================================================================================
/**
*
* @brief Model that holds and manages raw fiff data.
*/
class ANSHAREDSHARED_EXPORT FiffRawModel : public AbstractModel
{
    Q_OBJECT

public:
    typedef QSharedPointer<FiffRawModel> SPtr;            /**< Shared pointer type for FiffRawModel. */
    typedef QSharedPointer<const FiffRawModel> ConstSPtr; /**< Const shared pointer type for FiffRawModel. */

    //=========================================================================================================
    /**
    * Constructs a FiffRawModel object.
    */
    FiffRawModel(QObject *pParent = nullptr);

    //=========================================================================================================
    /**
    * Constructs a FiffRawModel object.
    */
    FiffRawModel(QFile& inFile,
                 qint32 iSamplesPerBlock,
                 qint32 iWindowSize,
                 qint32 iPaddingSize,
                 QObject *pParent = nullptr);

    //=========================================================================================================
    /**
    * Destructs a FiffRawModel.
    */
    ~FiffRawModel();

    //=========================================================================================================
    /**
    * Resets the model and reads new data from inFile
    *
    * @param inFile The file to read data from.
    */
    void initFiffData(QFile& inFile);

    //=========================================================================================================
    /**
    * Returns the data stored under the given role for the index.
    * Currently only Qt::DisplayRole is supported
    *
    * @param[in] index   The index that referres to the requested item.
    * @param[in] role    The requested role.
    */
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;

    //=========================================================================================================
    /**
    * Returns the item flags for the given index.
    *
    * @param[in] index   The index that referres to the requested item.
    */
    Qt::ItemFlags flags(const QModelIndex &index) const override;

    //=========================================================================================================
    /**
    * Returns the index for the item in the model specified by the given row, column and parent index.
    *
    * @param[in] row      The specified row.
    * @param[in] column   The specified column.
    * @param[in] parent   The parent index.
    */
    QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const override;

    //=========================================================================================================
    /**
    * Returns the parent index of the given index.
    * In this Model the parent index in always QModelIndex().
    *
    * @param[in] index   The index that referres to the child.
    */
    QModelIndex parent(const QModelIndex &index) const override;

    //=========================================================================================================
    /**
    * Returns the number of childeren for the parent node.
    *
    * @param[in] parent     The parent index.
    */
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;

    //=========================================================================================================
    /**
    * Returns the number of objects stored in the node.
    *
    * @param[in] parent     The index of the requested node.
    */
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;

    //=========================================================================================================
    /**
    * Returns true if parent has any children; otherwise returns false.
    *
    * @param[in] parent     The index of the parent node.
    */
    bool hasChildren(const QModelIndex &parent = QModelIndex()) const override;

    //=========================================================================================================
    /**
    * The type of this model (QEntityListModel)
    *
    * @return The type of this model (QEntityListModel)
    */
    inline MODEL_TYPE getType() const override;

    //=========================================================================================================
    /**
    * Return the first sample of the loaded Fiff file
    *
    * @return The first sample of the loaded Fiff file
    */
    inline qint32 firstSample() const;

    //=========================================================================================================
    /**
    * Returns the last sample of the loaded Fiff file
    *
    * @return The last sample of the loaded Fiff file
    */
    inline qint32 lastSample() const;

signals:

    //=========================================================================================================
    /**
    * This is emitted in order to cheat the constness of the ::data method
    *
    * @param[in] iCursorRequested Cursor that points to the requested sample
    */
    void startToLoadBlocks(qint32 iCursorRequested) const;

private slots:

    //=========================================================================================================
    /**
    * This method determines whether we need to load earlier or later blocks and then calls the suitable
    * method in the background.
    *
    * @param iCursorRequested Cursor that points to the requested sample
    */
    void onStartToLoadBlocks(qint32 iCursorRequested);

private:

    //=========================================================================================================
    /**
    * This is run concurrently
    *
    * @param[in] iCursorRequested Cursor that points to the requested sample
    */
    int loadEarlierBlocks(qint32 iCursorRequested);

    //=========================================================================================================
    /**
    * This is run concurrently
    *
    * @param[in] iCursorRequested Cursor that points to the requested sample
    */
    int loadLaterBlocks(qint32 iCursorRequested);

    //=========================================================================================================
    /**
    * This is run by the FutureWatcher when its finished
    *
    * @param[in] result Code value for the result
    */
    void postBlockLoad(int result);

private:

    QList<QPair<MatrixXd, MatrixXd>> m_lData;    /**< Data */

    qint32 m_iSamplesPerBlock;  /**< Number of samples per block */
    qint32 m_iWindowSize;       /**< Number of blocks per window */
    qint32 m_iPaddingSize;      /**< Number of blocks that are padded left and right */

    // this always points to the very first sample that is currently held (in the earliest block)
    qint32 m_iFiffCursorBegin;


    // concurrent reloading
    QFutureWatcher<int> m_blockLoadFutureWatcher;    /**< QFutureWatcher for watching process of reloading fiff data. */


    QSharedPointer<FIFFLIB::FiffIO> m_pFiffIO;
};


//*************************************************************************************************************
//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================

inline MODEL_TYPE FiffRawModel::getType() const
{
    return MODEL_TYPE::ANSHAREDLIB_FIFFRAW_MODEL;
}


//*************************************************************************************************************

inline qint32 FiffRawModel::firstSample() const {
    if(m_pFiffIO->m_qlistRaw.empty() == false)
        return m_pFiffIO->m_qlistRaw[0]->first_samp;
    else
    {
        qDebug() << "[FiffRawModel::firstSample] Raw list is empty, returning -1";
    }
}


//*************************************************************************************************************

inline qint32 FiffRawModel::lastSample() const {
    if(m_pFiffIO->m_qlistRaw.empty() == false)
        return m_pFiffIO->m_qlistRaw[0]->last_samp;
    else
    {
        qDebug() << "[FiffRawModel::lastSample] Raw list is empty, returning -1";
    }
}

} // namespace ANSHAREDLIB

#endif // ANSHAREDLIB_FIFFRAWMODEL_H
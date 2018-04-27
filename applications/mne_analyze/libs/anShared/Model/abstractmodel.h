//=============================================================================================================
/**
* @file     abstractmodel.h
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
* @brief     AbstractModel class declaration.
*
*/

#ifndef ANSHAREDLIB_ABSTRACTMODEL_H
#define ANSHAREDLIB_ABSTRACTMODEL_H


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../anshared_global.h"

//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>
#include <QAbstractItemModel>


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
* Description of what this class is intended to do (in detail).
*
* @brief Brief description of this class.
*/
class ANSHAREDSHARED_EXPORT AbstractModel : public QAbstractItemModel
{
    Q_OBJECT

public:
    typedef QSharedPointer<AbstractModel> SPtr;            /**< Shared pointer type for AbstractModel. */
    typedef QSharedPointer<const AbstractModel> ConstSPtr; /**< Const shared pointer type for AbstractModel. */

    //=========================================================================================================
    /**
    * Constructs a AbstractModel object. Simply pass potential parent object to super class.
    */
    AbstractModel(QObject *pParent = nullptr)
        : QAbstractItemModel(pParent) {}

    //=========================================================================================================
    /**
    * Default destructor.
    */
    virtual ~AbstractModel() = default;

    //=========================================================================================================
    /**
    * @brief The MODEL_TYPE enum lists all available model types.
    *        Naming convention: NAMESPACE_CLASSNAME_MODEL
    */
    enum MODEL_TYPE
    {
        FSLIB_SURFACE_MODEL
    };

    //=========================================================================================================
    /**
    * @brief getType Inherited by AbstractModel
    * @return The type of the respective subclasses
    */
    virtual inline MODEL_TYPE getType() const = 0;

    //=========================================================================================================
    // Inherited by QAbstractItemModel:
    virtual QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override = 0;
    virtual Qt::ItemFlags flags(const QModelIndex &index) const override = 0;
    virtual QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const override = 0;
    virtual QModelIndex parent(const QModelIndex &index) const override = 0;
    virtual int rowCount(const QModelIndex &parent = QModelIndex()) const override = 0;
    virtual int columnCount(const QModelIndex &parent = QModelIndex()) const override = 0;

protected:

private:

};


//*************************************************************************************************************
//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================


} // namespace ANSHAREDLIB

#endif // ANSHAREDLIB_ABSTRACTMODEL_H
//=============================================================================================================
/**
* @file     abstracttreeitem.h
* @author   Lorenz Esch <Lorenz.Esch@tu-ilmenau.de>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     November, 2015
*
* @section  LICENSE
*
* Copyright (C) 2015, Lorenz Esch and Matti Hamalainen. All rights reserved.
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
* @brief     AbstractTreeItem class declaration.
*
*/

#ifndef ABSTRACTTREEITEM_H
#define ABSTRACTTREEITEM_H

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "types.h"


//*************************************************************************************************************
//=============================================================================================================
// Qt INCLUDES
//=============================================================================================================

#include <QStandardItem>


//*************************************************************************************************************
//=============================================================================================================
// Eigen INCLUDES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// DEFINE NAMESPACE DISP3DLIB
//=============================================================================================================

namespace DISP3DLIB
{


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================


//=============================================================================================================
/**
* AbstractTreeItem provides AbstractTreeItem provides the basic tree item. This item should be used as a base class for all tree items throughout the disp3D library.
*
* @brief Provides the basic tree item.
*/
class AbstractTreeItem : public QObject, public QStandardItem
{
    Q_OBJECT;

public :
    //=========================================================================================================
    /**
    * Default constructor.
    *
    * @param[in] iType      The type of the item. See types.h for declaration and definition.
    * @param[in] text       The text of this item. This is also by default the displayed name of the item in a view.
    */
    AbstractTreeItem(const int& iType, const QString& text = "");
    virtual ~AbstractTreeItem();

    //=========================================================================================================
    /**
    * QStandardItem functions
    */
    QVariant data(int role = Qt::UserRole + 1) const;
    void setData(const QVariant& value, int role = Qt::UserRole + 1);
    int type() const;

    //=========================================================================================================
    /**
    * Returns all children of this item based on their type.
    *
    * @param[in] type    The type of the child items which should be looked for.
    *
    * @return           List with all found items.
    */
    QList<QStandardItem*> findChildren(const int& type);

    //=========================================================================================================
    /**
    * Returns all children of this item based on their text.
    *
    * @param[in] text    The text of the child items which should be looked for.
    *
    * @return           List with all found items.
    */
    QList<QStandardItem*> findChildren(const QString& text);

    //=========================================================================================================
    /**
    * Overloaded stream operator to add a child to this item based on a pointer.
    *
    * @param[in] newItem    The new item as a pointer.
    */
    AbstractTreeItem &operator<<(AbstractTreeItem* newItem);

    //=========================================================================================================
    /**
    * Overloaded stream operator to add a child to this item based on a reference.
    *
    * @param[in] newItem    The new item as a reference.
    */
    AbstractTreeItem &operator<<(AbstractTreeItem& newItem);

protected:
    int     m_iType;        /**< This item's type. */

signals:
    //=========================================================================================================
    /**
    * Emit this signal whenever this item's check state changed.
    *
    * @param[in] checkState     The current check state.
    */
    void checkStateChanged(const Qt::CheckState& checkState);
};

} //NAMESPACE DISP3DLIB

#endif // ABSTRACTTREEITEM_H
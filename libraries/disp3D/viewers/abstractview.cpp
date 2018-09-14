//=============================================================================================================
/**
* @file     abstractview.cpp
* @author   Lorenz Esch <Lorenz.Esch@tu-ilmenau.de>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     March, 2017
*
* @section  LICENSE
*
* Copyright (C) 2017, Lorenz Esch and Matti Hamalainen. All rights reserved.
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
* @brief    AbstractView class definition.
*
*/

//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "abstractview.h"

#include "../engine/view/view3D.h"
#include "../engine/model/data3Dtreemodel.h"
#include "../engine/control/control3dwidget.h"

#include <inverse/dipoleFit/dipole_fit_settings.h>
#include <inverse/dipoleFit/ecd_set.h>
#include <mne/mne_bem.h>
#include <disp/viewers/quickcontrolview.h>


//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QGridLayout>
#include <QSizePolicy>


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace DISP3DLIB;
using namespace DISPLIB;
using namespace INVERSELIB;
using namespace MNELIB;


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

AbstractView::AbstractView(QWidget* parent,
                           Qt::WindowFlags f)
: QWidget(parent, f)
, m_p3DView(View3D::SPtr(new View3D()))
, m_pData3DModel(Data3DTreeModel::SPtr(new Data3DTreeModel()))
{
    //Init 3D View
    m_p3DView->setModel(m_pData3DModel);
    m_p3DView->setFlag(Qt::FramelessWindowHint, true);

    QStringList slControlFlags;
    slControlFlags << "Data" << "View" << "Light";
    m_pControl3DView = Control3DWidget::SPtr(new Control3DWidget(this, slControlFlags));

    m_pControl3DView->init(m_pData3DModel, m_p3DView);

    createGUI();
}


//*************************************************************************************************************

AbstractView::~AbstractView()
{
}


//*************************************************************************************************************

QSharedPointer<DISP3DLIB::View3D> AbstractView::getView()
{
    return m_p3DView;
}


//*************************************************************************************************************

QSharedPointer<DISP3DLIB::Control3DWidget> AbstractView::getControlView()
{
    return m_pControl3DView;
}


//*************************************************************************************************************

QSharedPointer<DISP3DLIB::Data3DTreeModel> AbstractView::getTreeModel()
{
    return m_pData3DModel;
}


//*************************************************************************************************************

QPointer<QuickControlView> AbstractView::getQuickControl()
{
    return m_pQuickControlView;
}


//*************************************************************************************************************

void AbstractView::setQuickControlWidgets(const QList<QWidget*>&lControlWidgets)
{
    if(m_pQuickControlView) {
        for(int i = 0; i < lControlWidgets.size(); i++) {
            m_pQuickControlView->addGroupBox(lControlWidgets.at(i), lControlWidgets.at(i)->windowTitle());
        }
    }
}


//*************************************************************************************************************

void AbstractView::createGUI()
{
    m_pQuickControlView = new QuickControlView("Network View", Qt::Widget, this, false);
    m_pQuickControlView->setVisiblityHideOpacityClose(false);

    //Create widget GUI
    m_pQuickControlView->addGroupBox(m_pControl3DView.data(), "3D View");

    QWidget *pWidgetContainer = QWidget::createWindowContainer(m_p3DView.data(), this, Qt::Widget);
    pWidgetContainer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    pWidgetContainer->setMinimumSize(400,400);

    QGridLayout* pMainLayoutView = new QGridLayout();
    pMainLayoutView->addWidget(pWidgetContainer,0,0);
    pMainLayoutView->addWidget(m_pQuickControlView.data(),0,1);


    this->setLayout(pMainLayoutView);
}

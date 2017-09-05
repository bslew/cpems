/********************************************************************************
** Form generated from reading ui file 'progressbar.ui'
**
** Created: Fri May 14 23:32:34 2010
**      by: Qt User Interface Compiler version 4.5.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_PROGRESSBAR_H
#define UI_PROGRESSBAR_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QProgressBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_progressBar
{
public:
    QGridLayout *gridLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *node;
    QProgressBar *bar;
    QLabel *range;

    void setupUi(QWidget *progressBar)
    {
        if (progressBar->objectName().isEmpty())
            progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->resize(216, 29);
        gridLayout = new QGridLayout(progressBar);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setSizeConstraint(QLayout::SetMinimumSize);
        gridLayout->setContentsMargins(-1, 0, -1, 0);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        node = new QLabel(progressBar);
        node->setObjectName(QString::fromUtf8("node"));

        horizontalLayout->addWidget(node);

        bar = new QProgressBar(progressBar);
        bar->setObjectName(QString::fromUtf8("bar"));
        bar->setValue(24);

        horizontalLayout->addWidget(bar);

        range = new QLabel(progressBar);
        range->setObjectName(QString::fromUtf8("range"));

        horizontalLayout->addWidget(range);


        gridLayout->addLayout(horizontalLayout, 0, 0, 1, 1);


        retranslateUi(progressBar);

        QMetaObject::connectSlotsByName(progressBar);
    } // setupUi

    void retranslateUi(QWidget *progressBar)
    {
        progressBar->setWindowTitle(QApplication::translate("progressBar", "Form", 0, QApplication::UnicodeUTF8));
        node->setText(QApplication::translate("progressBar", "node:", 0, QApplication::UnicodeUTF8));
        range->setText(QApplication::translate("progressBar", "from/to", 0, QApplication::UnicodeUTF8));
        Q_UNUSED(progressBar);
    } // retranslateUi

};

namespace Ui {
    class progressBar: public Ui_progressBar {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PROGRESSBAR_H
